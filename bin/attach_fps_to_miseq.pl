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
use LIMS2::Model::Util::Miseq qw(wells_generator);
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
    my ($lims2_model, $row, $fp_well) = @_;

    say 'Creating process - ' . $row->{'Cell Number'} . ':' . $fp_well->name . ' -> ' 
        . $row->{'Miseq Plate'} . ':' . $fp_well->name;

    my $adjustment = $row->{'Index Min'};
    my $end = $row->{'Index Max'};

    my $well_name = $fp_well->name;
    if ($adjustment) {
        my $int_well = wells_generator(1);

        my $well_pos = $int_well->{ $fp_well->name };

        my $start_pos = $adjustment % 96;


        my $end_pos = $end % 96;
        if ($end_pos == 0) {
            $end_pos = 96;
        }
        my $quadrant = int($adjustment / 96);

        if ($well_pos > $end_pos || $well_pos < $start_pos) {
            return;
        }

        my @name_well = wells_generator();

        my $miseq_adjust = ( $quadrant * 96 ) + $well_pos;
        $well_name = $name_well[$miseq_adjust];
    }

    my $miseq_well_rs = $lims2_model->schema->resultset('Well')->find({
        name => $well_name,
        'plate.name' => $row->{'Miseq Plate'},
    }, { 
        prefetch => 'plate'
    });

    unless ($miseq_well_rs) {
        return {
            fp          => $row->{'FP Plate'},
            well        => $fp_well->name,
            miseq_well  => $well_name,
            miseq_plate => $row->{'Miseq Plate'},
            reason      => 'No Miseq Child well found',
        };
    }

    my @input = map { $_->as_hash } $fp_well;
    my @output = ($miseq_well_rs->as_hash);

    my $process = {
        type            => 'miseq_no_template',
        input_wells     => \@input,
        output_wells    => \@output,
    };

    $lims2_model->create_process($process);


    my $ph_name = $row->{'Miseq Plate'} . '_FP';
    print "\nLooking for placeholder - $ph_name";
    my $placeholder_fp = $lims2_model->schema->resultset('Plate')->find({ name => $ph_name });

    if ($placeholder_fp) {
        print "\nChecking for Placeholder well";
        my $placeholder_well = $lims2_model->schema->resultset('Well')->find({
            name => $well_name,
            plate_id => $placeholder_fp->id
        });
        if ($placeholder_well) {
            my @ph_processes = map { $_->id } $placeholder_well->child_processes;
            foreach my $ph_process (@ph_processes) {
                print "\nDeleting place-holder process: $ph_process";
                $lims2_model->delete_process({ id => $ph_process });
            }
            say "\nDeleting place-holder well: $placeholder_fp " . $placeholder_well->name;
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
    if ($row->{'Miseq Plate'} eq 'No Miseq') {
        next;
    }

    my $fp_plate = $lims2_model->schema->resultset('Plate')->find({ name => $row->{'FP Plate'} });
    my $miseq_exp = $lims2_model->schema->resultset('MiseqExperiment')->search({ name => { like => $row->{'Cell Number'} . '%' } });
    while (my $exp = $miseq_exp->next) {
        unless ( $exp->parent_plate_id && $exp->experiment_id ) {
            say "Updating Miseq Exp: " . $exp->id . ' Parent: ' . $row->{'FP Plate'};
            if ($fp_plate) {
                $lims2_model->update_miseq_experiment({
                    id                  => $exp->id,
                    parent_plate_id     => $fp_plate->id,
                    experiment_id       => $row->{'Experiment ID'},
                });
            } else {
                push @errors, "Cant find plate: " . $row->{'FP Plate'};
            }
        }
    }

    unless ($fp_plate) {
        push (@errors, "No FP plate found - " . $row->{'Cell Number'});
        next;
    }

    my @fp_wells = map { $_ } $fp_plate->wells;
    foreach my $well (@fp_wells) {
        my @child_wells = map { $_->well } map { $_->process->process_output_wells } $well->process_input_wells->all;

        my $error;
        unless (grep { $_->plate->name eq $row->{'Miseq Plate'} && $_->name eq $well->name } @child_wells) {
            print Dumper $row->{'Miseq Plate'} . $well->name;
            $error = create_relation($lims2_model, $row, $well);
        }
        if ($error) {
            push @errors, $error;
        }
    }
}

my $json = JSON->new;
my $error_json = encode_json \@errors;
open(my $mfh, '>', 'fp_error.json');
print $mfh Dumper $error_json;
close $mfh;
