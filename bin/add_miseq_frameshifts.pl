#! /usr/bin/perl

use strict;
use warnings;

use LIMS2::Model;

my $model = LIMS2::Model->new({ user => 'lims2' });

sub get_plate_name {
    my $num = shift;
    my $plate_name;
    if ($num < 100) {
        $plate_name = "Miseq_0$num";
    } else {
        $plate_name = "Miseq_$num";
    }
    return $plate_name;
}

sub get_total_reads {
    my @miseq_alleles_frequencies = @_;
    my $total_reads = 0;
    foreach my $miseq_alleles_frequency (@miseq_alleles_frequencies) {
        $total_reads += $miseq_alleles_frequency->n_reads;
    }
    return $total_reads
}

sub check_read {
    my $read = shift;
    if ($read->nhej == 1) {
        if (($read->n_deleted + $read->n_inserted) % 3) {
            return 1;
        }
    }
    return;
}

sub check_frameshift {
    my @miseq_alleles_frequencies = @_;
    my $total_reads = get_total_reads(@miseq_alleles_frequencies);
    if (($miseq_alleles_frequencies[2]->n_reads / $total_reads) < 0.05) {
        my $read_1_fs = check_read($miseq_alleles_frequencies[0]);
        my $read_2_fs = check_read($miseq_alleles_frequencies[1]);
        if ($read_1_fs or $read_2_fs) {
            return 1;
        }
    }
    return;
}

foreach my $num (82..112) {
    my $plate_name = get_plate_name($num);
    my $plate_id = $model->schema->resultset('Plate')->find({name => $plate_name})->id;
    my $miseq_plate = $model->schema->resultset('MiseqPlate')->find({plate_id => $plate_id});
    foreach my $miseq_experiment ($miseq_plate->miseq_experiments) {
        foreach my $miseq_well_experiment ($miseq_experiment->miseq_well_experiments) {
            my @miseq_alleles_frequencies = $model->schema->resultset('MiseqAllelesFrequency')->search(
                {miseq_well_experiment_id => $miseq_well_experiment->id},
                {order_by => {'-desc' => 'n_reads'}}
            );
            if (scalar @miseq_alleles_frequencies >= 3) {
                if (check_frameshift(@miseq_alleles_frequencies)) {
                    $model->update_miseq_well_experiment({
                            id => $miseq_well_experiment->id,
                            frameshifted => 1});
                }
            }
        }
    }
}
