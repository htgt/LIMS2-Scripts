#! /usr/bin/perl

use strict;
use warnings;

use LIMS2::Model;
use Try::Tiny;
use Getopt::Long;
use Data::Dumper;
use Moose;
use LIMS2::REST::Client;
use Bio::Perl qw( revcom );
use Text::CSV;
use POSIX qw(strftime);
use feature qw(say);
use List::MoreUtils qw(uniq);
use JSON;

sub extract_rows {
    my ($csv, $fh, @data) = @_;

    my $headers = $csv->getline($fh);
    $csv->column_names( @{ $headers} );

    while (my $row = $csv->getline_hr($fh)) {
        push (@data, $row);
    }

    return @data;
}

sub create_relation {
    my ($lims2_model, $row, $miseq_well) = @_;

    say 'Creating process - ' . $row->{PIQ_Plate} . ':' . $row->{PIQ_Well} . ' -> ' 
        . $row->{Miseq_Exp_Plate} . ':' . $miseq_well->name . '~' . $row->{Miseq_Experiment};

    my $input_rs = $lims2_model->schema->resultset('Well')->find({
        name => $row->{PIQ_Well},
        'plate.name' => $row->{PIQ_Plate}
    }, { 
        prefetch => 'plate'
    });
    my @input;
    if ($input_rs) {
        @input = map { $_->as_hash } $input_rs;
    } else {
        return {
            piq         => $row->{PIQ_Plate},
            piq_well    => $row->{PIQ_Well},
            miseq_plate => $row->{Miseq_Exp_Plate},
            miseq_well  => $miseq_well->name,
            miseq_exp   => $row->{Miseq_Experiment},
            reason      => 'No parent well found',
        };
    }
    my @output = ($miseq_well->as_hash);

    my $process = {
        type            => 'miseq_no_template',
        input_wells     => \@input,
        output_wells    => \@output,
    };

    $lims2_model->create_process($process);

    my $ph_name = $row->{Miseq_Exp_Plate} . '_FP';
    print "\nLooking for placeholder - $ph_name";
    my $placeholder_fp = $lims2_model->schema->resultset('Plate')->find({ name => $ph_name });

    if ($placeholder_fp) {
        print "\nChecking for Placeholder well";
        my $placeholder_well = $lims2_model->schema->resultset('Well')->find({
            name => $miseq_well->name,
            plate_id => $placeholder_fp->id
        });
        if ($placeholder_well) {
            my @ph_processes = map { $_->id } $placeholder_well->child_processes;
            foreach my $ph_process (@ph_processes) {
                print "\nDeleting place-holder process: $ph_process";
                $lims2_model->delete_process({ id => $ph_process });
            }
            print "\nDeleting place-holder well: $placeholder_fp $process->{output_wells}[0]->{well_name}";
            $lims2_model->delete_well({ id => $placeholder_well->id });
        }
    }


    return;
}


my ($file, $qc_bool);

GetOptions(
    'file=s'  => \$file,
)
or die usage_message();;

my $lims2_model = LIMS2::Model->new( user => 'lims2' );
my $csv = Text::CSV->new({ binary => 1 }) or die "Cannot use CSV: " . Text::CSV->error_diag();
my @data;

my $fh;
open $fh, "<:encoding(utf8)", $file or die;
@data = extract_rows($csv, $fh, @data);
close $fh;

my @errors;
foreach my $row (@data) {

    my $miseq_well_exps = $lims2_model->schema->resultset('MiseqWellExperiment')->search({ miseq_exp_id => $row->{Miseq_Exp_ID} },{ prefetch => 'well'});
    my $piq_plate = $lims2_model->schema->resultset('Plate')->find({ name => $row->{PIQ_Plate} });
    my $piq_well = $lims2_model->schema->resultset('Well')->find({ 
        name => $row->{PIQ_Well},
        'plate.name' => $row->{PIQ_Plate}
    },{
        prefetch => 'plate'
    });

    if ($piq_well) {
        my @exps = map { $_->id } $piq_well->experiments_pipelineII;

        if (scalar @exps > 1) {
            push @errors, {
                piq         => $row->{PIQ_Plate},
                piq_well    => $row->{PIQ_Well},
                miseq_plate => $row->{Miseq_Exp_Plate},
                miseq_exp   => $row->{Miseq_Experiment},
                reason      => 'Multiple experiments inherited',
                exps        => \@exps,
            };;
        }
        unless ( $miseq_well_exps->first->miseq_exp->parent_plate_id && $miseq_well_exps->first->miseq_exp->experiment_id ) {
            $lims2_model->update_miseq_experiment({
                id                  => $miseq_well_exps->first->miseq_exp->id,
                parent_plate_id     => $piq_plate->id,
                experiment_id       => $exps[0],
            })
        }
    }
    print Dumper $miseq_well_exps->first->miseq_exp->as_hash;


    while (my $well_exp = $miseq_well_exps->next) {
        my $well = $well_exp->well;
        say $row->{PIQ_Plate} . ':' . $row->{PIQ_Well} . ' -> ' 
        . $row->{Miseq_Exp_Plate} . ':' . $well->name . ' ~ ' . $row->{Miseq_Experiment};

        my @parent_wells = map { $_->well } map { $_->process->process_input_wells } $well->process_output_wells->all;

        my $error;
        unless (grep { $_->plate->name eq $row->{PIQ_Plate} && $_->name eq $row->{PIQ_Well} } @parent_wells) {
            $error = create_relation($lims2_model, $row, $well);
        }
        if ($error) {
            push @errors, $error;
        }
    }
}

my $json = JSON->new;
my $error_json = encode_json \@errors;
open(my $mfh, '>', 'piq_error.json');
print $mfh Dumper $error_json;
close $mfh;