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

sub extract_rows {
    my ($csv, $fh, @data) = @_;

    my $headers = $csv->getline($fh);
    $csv->column_names( @{ $headers} );

    while (my $row = $csv->getline_hr($fh)) {
        push (@data, $row);
    }

    return @data;
}

sub map_well_names {
    my ($prefix, @nums) = @_;

    return map { sprintf("%s%02d", $prefix, $_) } @nums;
}

my ($file, $qc_bool);

GetOptions(
    'file=s'  => \$file,
    'primary_qc' => \$qc_bool,
)
or die usage_message();;

my $lims2_model = LIMS2::Model->new( user => 'lims2' );
my $csv = Text::CSV->new({ binary => 1 }) or die "Cannot use CSV: " . Text::CSV->error_diag();
my @data;

my $fh;
open $fh, "<:encoding(utf8)", $file or die;
@data = extract_rows($csv, $fh, @data);
close $fh;

foreach my $relation (@data) {
    say "Importing $relation->{miseq_plate_name}_$relation->{miseq_experiment}";

    my $parent_plate = $lims2_model->schema->resultset('Plate')->find({ name => $relation->{parent_plate_name} });
    my $miseq_plate = $lims2_model->schema->resultset('Plate')->find({ name => $relation->{miseq_plate_name} });

    my $placeholder_fp = $lims2_model->schema->resultset('Plate')->find({ name => $relation->{miseq_plate_name} . "_FP" });

    unless ($parent_plate && $miseq_plate) {
        next;
    }

    say "Updating miseq experiment: $relation->{miseq_experiment}";
    my $miseq_details = $lims2_model->schema->resultset('MiseqPlate')->find({ plate_id => $miseq_plate->id });
    my $miseq_experiment = $lims2_model->schema->resultset('MiseqExperiment')->find({
        miseq_id => $miseq_details->id,
        name => $relation->{miseq_experiment},
    });

    $lims2_model->update_miseq_experiment({ 
        id              => $miseq_experiment->id,
        experiment_id   => $relation->{experiment_id},
        parent_plate_id => $parent_plate->id,
        gene            => $relation->{trivial_id},
    });

    print Dumper $miseq_experiment = $lims2_model->schema->resultset('MiseqExperiment')->find({
            miseq_id => $miseq_details->id,
            name => $relation->{miseq_experiment},
    })->as_hash;
    
    my @wells;
    my $well_str = $relation->{wells};

    if ($qc_bool) {
        my ($start_letter, $start_num, $end_letter, $end_num) = ($well_str =~ /\d*\ ?\(([A-P])([0-9]{2})\-([A-P])([0-9]{2})\)/g); #Capture start and end wells
        my @letters = $start_letter .. $end_letter;
        my @numbers = $start_num .. $end_num;
        @wells = map { map_well_names($_, @numbers) } @letters;
    }

    say "Connecting wells";
    foreach my $well (@wells) {
        my $parent_well = $lims2_model->schema->resultset('Well')->find({ name => $well, plate_id => $parent_plate->id });
        my $miseq_well = $lims2_model->schema->resultset('Well')->find({ name => $well, plate_id => $miseq_plate->id });

        my %input_child_wells = map { $_->id => 1 } $parent_well->child_wells;

        if (!($input_child_wells{$miseq_well->id})) {
            my $process = {
                type => 'miseq_no_template',
                input_wells => [{
                    plate_name  => $parent_plate->name,
                    well_name   => $well,
                }],
                output_wells => [{
                    plate_name  => $miseq_plate->name,
                    well_name   => $well,               
                }],
            };
            my $process_rs = $lims2_model->create_process($process);
        } else {
            say "Relation already exists. $parent_well -> $miseq_well";
        }
        if ($placeholder_fp) {
            my $placeholder_well = $lims2_model->schema->resultset('Well')->find({ name => $well, plate_id => $placeholder_fp->id });
            if ($placeholder_well) {
                my @ph_processes = map { $_->id } $placeholder_well->child_processes;
                foreach my $ph_process (@ph_processes) {
                    say "Deleting place-holder process: $ph_process";
                    $lims2_model->delete_process({ id => $ph_process });
                }
                say "Deleting place-holder well: $placeholder_fp $well";
                $lims2_model->delete_well({ id => $placeholder_well->id });
            }
        }
    }

    my @placeholder_wells = map { $_->id } $placeholder_fp->wells;
    unless (@placeholder_wells) {
        say "Deleting placeholder plate: $placeholder_fp->name";
        $lims2_model->delete_plate({ id => $placeholder_fp->id });
    }
}


1;
