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

sub clean_up_placeholders {
    my ($lims2_model, $placeholder_fp, $well) = @_;
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

    my @placeholder_wells = map { $_->id } $placeholder_fp->wells;
    unless (@placeholder_wells) {
        say "Deleting placeholder plate: $placeholder_fp->name";
        $lims2_model->delete_plate({ id => $placeholder_fp->id });
    }
}

sub clean_up_tt_ph {
    my ($lims2_model, $tt_num, $well, $mw_id) = @_;

    my $tt_plate_name = "MiSeq_TT_FP_$tt_num";
    my $tt_well_rs = $lims2_model->schema->resultset('Well')->find({ 
        name => $well,
        'plate.name' => $tt_plate_name,
    },
    {
        join => 'plate'
    });

    if ($tt_well_rs) {
        #my @ph_processes = map { $_->id } $placeholder_well->child_processes;
        my $relation = $lims2_model->schema->resultset('Process')->find({
            'process_input_wells.well_id' => $tt_well_rs->id,
            'process_output_wells.well_id' => $mw_id,
        }, {
            join => ['process_output_wells', 'process_input_wells'] 
        });
        if ($relation) {
            say "Deleting place-holder process: $relation";
            $lims2_model->delete_process({ id => $relation->id });
        }
        my $tt_plate = $tt_well_rs->plate;
        my @tt_pro = map { $_->id } $tt_well_rs->child_processes;
        unless (@tt_pro) {
            say "Deleting place-holder well: $tt_well_rs";
            $lims2_model->delete_well({ id => $tt_well_rs->id });
        }
        my @placeholder_wells = map { $_->id } $tt_plate->wells;
        unless (@placeholder_wells) {
            say "Deleting placeholder plate: $tt_plate->name";
            $lims2_model->delete_plate({ id => $tt_plate->id });
        }
    }
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

my $plates;
my $row = 1;
my @dis_fp;
my @dis_piq;
foreach my $relation (@data) {
    my $primary_key = "$relation->{HUPFP_Plate}_$relation->{Well_Range}_$relation->{'1st_Miseq_Plate'}";
    my $primary_hash = {
        parent_plate        => $relation->{HUPFP_Plate},
        miseq_plate         => $relation->{'1st_Miseq_Plate'},
        miseq_experiment    => $relation->{'1st_Miseq_Experiment'},
        well_range          => $relation->{Well_Range},
        exp_id              => $relation->{Experiment_ID},
        trivial_id          => $relation->{Experiment},
    };

    if ( my $curr = $plates->{FP}->{$primary_key} ) {
        foreach my $key (keys %{ $primary_hash }) {
            if ($curr->{$key} ne $primary_hash->{$key}) {
                say "Row $row: $primary_key 1stQC - $key data discrepancy { Found: $curr->{$key}, Row: $primary_hash->{$key} }";
                push (@dis_fp, $primary_key);
            }
        }
    } else {
        $plates->{FP}->{$primary_key} = $primary_hash;
    }

    if ($relation->{'2nd_Miseq_Plate'} ) {
        my $secondary_key = "$relation->{HUEDQ_Plate}_$relation->{PIQ_Well}_$relation->{'2nd_Miseq_Plate'}";
        my $secondary_hash = {
            parent_plate        => $relation->{HUEDQ_Plate},
            parent_well         => $relation->{PIQ_Well},
            miseq_plate         => $relation->{'2nd_Miseq_Plate'},
            miseq_well          => $relation->{Miseq_Well},
            miseq_experiment    => $relation->{'2nd_Miseq_Experiment'},
            exp_id              => $relation->{Experiment_ID},
            trivial_id          => $relation->{Experiment},
            primary_qc          => $relation->{'1st_Miseq_Plate'},
        };
        if ( my $curr = $plates->{PIQ}->{$secondary_key} ) {
            foreach my $key (keys %{ $secondary_hash }) {
                if ($curr->{$key} ne $secondary_hash->{$key}) {
                    say "Row $row: $secondary_key 2ndQC - $key data discrepancy { Found: $curr->{$key}, Row: $secondary_hash->{$key} }";
                    push (@dis_piq, $secondary_key);
                }
            }
        } else {
            $plates->{PIQ}->{$secondary_key} = $secondary_hash; 
        }
    }
    $row++;
}

@dis_fp = uniq @dis_fp;
@dis_piq = uniq @dis_piq;

map { delete $plates->{FP}->{$_} } @dis_fp;
map { delete $plates->{PIQ}->{$_} } @dis_piq;
my @skipped;

foreach my $fp_key (keys %{ $plates->{FP} }) {
    my $fp_dets = $plates->{FP}->{$fp_key};
    say "Importing $fp_dets->{miseq_plate}_$fp_dets->{miseq_experiment}";

    my $parent_plate = $lims2_model->schema->resultset('Plate')->find({ name => $fp_dets->{parent_plate} });
    my $miseq_plate = $lims2_model->schema->resultset('Plate')->find({ name => $fp_dets->{miseq_plate} });

    my $placeholder_fp = $lims2_model->schema->resultset('Plate')->find({ name => $fp_dets->{miseq_plate} . "_FP" });

    unless ($parent_plate && $miseq_plate) {
        my $err = "$fp_dets->{parent_plate}-$fp_dets->{miseq_plate} not found.";
        say $err;
        push(@skipped, $err);
        next;
    }

    say "Updating miseq experiment: $fp_dets->{miseq_experiment}";
    my $miseq_details = $lims2_model->schema->resultset('MiseqPlate')->find({ plate_id => $miseq_plate->id });
    my $miseq_experiment = $lims2_model->schema->resultset('MiseqExperiment')->find({
        miseq_id => $miseq_details->id,
        name => $fp_dets->{miseq_experiment},
    });
    $lims2_model->update_miseq_experiment({ 
        id              => $miseq_experiment->id,
        experiment_id   => $fp_dets->{exp_id},
        parent_plate_id => $parent_plate->id,
        gene            => $fp_dets->{trivial_id},
    });
    print Dumper $miseq_experiment = $lims2_model->schema->resultset('MiseqExperiment')->find({
            miseq_id => $miseq_details->id,
            name => $fp_dets->{miseq_experiment},
    })->as_hash;
    
    my @wells;
    my $well_str = $fp_dets->{well_range};

    my ($start_letter, $start_num, $end_letter, $end_num) = ($well_str =~ /([A-P])([0-9]{2})\-([A-P])([0-9]{2})/g); #Capture start and end wells
    my @letters = $start_letter .. $end_letter;
    my @numbers = $start_num .. $end_num;
    @wells = map { map_well_names($_, @numbers) } @letters;


    say "Connecting wells";
    foreach my $well (@wells) {
        my $parent_well = $lims2_model->schema->resultset('Well')->find({ name => $well, plate_id => $parent_plate->id });
        my $miseq_well = $lims2_model->schema->resultset('Well')->find({ name => $well, plate_id => $miseq_plate->id });

        unless($parent_well) {
            my $err = "Skipped $fp_dets->{parent_plate}:$well - No Parent Well record found";
            say $err;
            push(@skipped, $err);
            next;
        }

        unless($miseq_well) {
            my $err = "Skipped $fp_dets->{miseq_plate}:$well - No Miseq Well record found";
            say $err;
            push(@skipped, $err);
            next;
        }

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
            clean_up_placeholders($lims2_model, $placeholder_fp, $miseq_well->name);
        }
        my ($tt_num) = $fp_dets->{parent_plate} =~ /HUPFPPlate(\d+)\w?/;
        if ($tt_num) {
            clean_up_tt_ph($lims2_model, $tt_num, $well, $miseq_well->id);
        }
    }

}

foreach my $piq_key (keys %{ $plates->{PIQ} }) {
    my $piq_dets = $plates->{PIQ}->{$piq_key};
    say "Importing $piq_dets->{miseq_plate}_$piq_dets->{miseq_experiment}";

    my $parent_plate = $lims2_model->schema->resultset('Plate')->find({ name => $piq_dets->{parent_plate} });
    my $miseq_plate = $lims2_model->schema->resultset('Plate')->find({ name => $piq_dets->{miseq_plate} });

    my $placeholder_fp = $lims2_model->schema->resultset('Plate')->find({ name => $piq_dets->{miseq_plate} . "_FP" });

    unless ($parent_plate && $miseq_plate) {
        say "$parent_plate-$miseq_plate not found.";
        next;
    }

    say "Updating miseq experiment: $piq_dets->{miseq_experiment}";
    my $miseq_details = $lims2_model->schema->resultset('MiseqPlate')->find({ plate_id => $miseq_plate->id });
    my $miseq_experiment = $lims2_model->schema->resultset('MiseqExperiment')->find({
        miseq_id => $miseq_details->id,
        name => $piq_dets->{miseq_experiment},
    });
    $lims2_model->update_miseq_experiment({ 
        id              => $miseq_experiment->id,
        experiment_id   => $piq_dets->{exp_id},
        parent_plate_id => $parent_plate->id,
        gene            => $piq_dets->{trivial_id},
    });
    print Dumper $miseq_experiment = $lims2_model->schema->resultset('MiseqExperiment')->find({
            miseq_id => $miseq_details->id,
            name => $piq_dets->{miseq_experiment},
    })->as_hash;
    
    say "Connecting well";
        
    my $parent_well = $lims2_model->schema->resultset('Well')->find({ name => $piq_dets->{parent_well}, plate_id => $parent_plate->id });
    my $miseq_well = $lims2_model->schema->resultset('Well')->find({ name => $piq_dets->{miseq_well}, plate_id => $miseq_plate->id });

    unless($parent_well) {
        my $err = "Skipped $piq_dets->{parent_plate}:$piq_dets->{parent_well} - No Parent Well record found";
        say $err;
        push(@skipped, $err);
        next;
    }

    unless($miseq_well) {
        my $err = "Skipped $piq_dets->{miseq_plate}:$piq_dets->{miseq_well} - No Miseq Well record found";
        say $err;
        push(@skipped, $err);
        next;
    }

    my %input_child_wells = map { $_->id => 1 } $parent_well->child_wells;

    if (!($input_child_wells{$miseq_well->id})) {
        my $process = {
            type => 'miseq_no_template',
            input_wells => [{
                plate_name  => $parent_plate->name,
                well_name   => $piq_dets->{parent_well},
            }],
            output_wells => [{
                plate_name  => $miseq_plate->name,
                well_name   => $piq_dets->{miseq_well},               
            }],
        };
        my $process_rs = $lims2_model->create_process($process);
    } else {
        say "Relation already exists. $parent_well -> $miseq_well";
    }

    if ($placeholder_fp) {
        clean_up_placeholders($lims2_model, $placeholder_fp, $miseq_well->name);
    }

}
say "Skipped values";
print Dumper @skipped;

1;
