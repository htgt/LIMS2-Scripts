#!/usr/bin/env perl

use strict;
use warnings;

use LIMS2::Model;
use LIMS2::Model::Util::Miseq qw(convert_index_to_well_name);

my $schema = LIMS2::Model->new( user => 'lims2' )->schema;
my $plate_rs = $schema->resultset('Plate');
my $miseq_plate_rs = $schema->resultset('MiseqPlate');
my $miseq_exp_rs = $schema->resultset('MiseqExperiment');
my $miseq_well_exp_rs = $schema->resultset('MiseqWellExperiment');
my $well_rs = $schema->resultset('Well');
my $well_not_found = 1;

sub get_all_miseq_plates {
    my @all_miseq_plates = map {$_->as_hash} $plate_rs->search({ type_id => 'MISEQ'});
    return @all_miseq_plates;
}

sub get_miseq_exps {
    my $plate_id = shift;
    #my $miseq_plate_id = $miseq_plate_rs->find({ plate_id => $plate_id })->id;
    my $miseq_plate_id = $miseq_plate_rs->find({ plate_id => $plate_id }, {column => 'id'})->id;
    my @miseq_exps = map {$_->as_hash} $miseq_exp_rs->search({ miseq_id => $miseq_plate_id });
    return @miseq_exps;
}

sub get_well_ids {
    my $exp_id = shift;
    #my @well_ids = map {$_->well_id} $miseq_well_exp_rs->search({ miseq_exp_id => $exp_id});
    my @well_ids = map {$_->well_id} $miseq_well_exp_rs->search({ miseq_exp_id => $exp_id}, {column => 'well_id'});
    return @well_ids;
}

sub get_well_name {
    my $well_id = shift;
    #my $well_name = $well_rs->find({ id => $well_id })->name;
    my $well_name = $well_rs->find({ id => $well_id }, {column => 'name'})->name;
    return $well_name || $well_not_found;
}

sub construct_file_path {
    my ($plate_name, $i, $exp_name) = @_;
    my $file_path = "/warehouse/team229_wh01/lims2_managed_miseq_data/$plate_name/S${i}_exp$exp_name/CRISPResso_on_${i}_S${i}_L001_R1_001_${i}_S${i}_L001_R2_001/Quantification_of_editing_frequency.txt";
    return $file_path;
}

sub find_mismatch {
    my ($file_path, $existing_well_names, $well_name_to_check) = @_;
    if (-e $file_path) {
        if (! exists($existing_well_names->{$well_name_to_check})) {
            return 1;
        }
    }
    return;
}

sub create_report_file {
    my ($plate_name, $mismatches) = @_;
    my $report_filename = "${plate_name}_lost_wells.txt";
    open(my $report_fh, '>', $report_filename) or die "Couldn't open $report_filename: $!";
    write_report($report_fh, $mismatches);
    close $report_fh;
    return;
}

sub write_report {
    my ($report_fh, $mismatches) = @_;
    while (my ($exp, $mismatch_wells) = each %$mismatches) {
        print $report_fh "$exp: @$mismatch_wells\n";
    }
    return;
}

#my @plates = get_all_miseq_plates();
my $plate = $plate_rs->find({'name' => 'Miseq_010'})->as_hash;
#foreach my $plate (@plates) {
    my %mismatches = ();
    print "Scanning $plate->{name}\n";
    my @exps = get_miseq_exps($plate->{id});
    foreach my $exp (@exps) {
        my @mismatch_wells = ();
        my @well_ids = get_well_ids($exp->{id});
        my %existing_well_names = map {get_well_name($_) => 1} @well_ids;
        foreach (my $i = 1; $i < 385; $i++) {
            my $file_path = construct_file_path($plate->{name}, $i, $exp->{name});
            my $well_name_to_check = convert_index_to_well_name($i);
            if (find_mismatch($file_path, \%existing_well_names, $well_name_to_check)) {
                push @mismatch_wells, $well_name_to_check;
            }
        }
        if (@mismatch_wells) {
            $mismatches{$exp->{name}} = \@mismatch_wells;
        }
    }
    if (keys %mismatches) {
        create_report_file($plate->{name}, \%mismatches);
    }
#}
