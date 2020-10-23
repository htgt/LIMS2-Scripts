#!/usr/bin/env perl

use strict;
use warnings;

use File::Slurp;
use Text::CSV;
use LIMS2::Model;

my $model = LIMS2::Model->new( user => 'lims2' );

my @output_data; 
my @headers = ('Label ID', 'Clone', 'CRISPR Sequence', 'WGE ID', 'Strand', 'Forward Primer', 'Reverse Primer', 'Miseq Plate');
push @output_data, \@headers;
my $file = $ARGV[0];
my @input_data = read_file($file, {chomp => 1});
foreach my $line (@input_data) {
    my @output_row = ();
    my ($fp_name, $well_name) = split ',', $line;
    ($fp_name) = $fp_name =~ m/([A-Za-z0-9_]+)/;
    ($well_name) = $well_name =~ m/([A-Z0-9]+)/;
    push @output_row, $fp_name, $well_name;
    my @blank_row = ($fp_name, $well_name, '', '', '', '', '');
    my $epd_name = get_parent_name($fp_name);
    if (! $epd_name) {
        print "Could not get epd name for $fp_name!\n";
        push @blank_row, '';
        push @output_data, \@blank_row;
        next;
    }
    my $miseq_plate_name = '';
    my $miseq_exp_id = get_miseq_exp_id($epd_name);
    if ($miseq_exp_id) {
        $miseq_plate_name = get_miseq_plate_name($miseq_exp_id);
    }
    push @blank_row, $miseq_plate_name;
    my $ep_name = get_parent_name($epd_name);
    if (! $ep_name) {
        print "Could not get ep name for $fp_name!\n";
        push @output_data, \@blank_row;
        next;
    }
    my $ep_id = get_ep_id($ep_name);
    if (! $ep_id) {
        print "Could not get ep id for $fp_name!\n";
        push @output_data, \@blank_row;
        next;
    }
    my $well_id = get_well_id($ep_id, $well_name);
    if (! $well_id) {
        print "Could not get well id for $fp_name, $well_name!\n";
        push @output_data, \@blank_row;
        next;
    }
    my $process_id = get_process_id($well_id);
    if (! $process_id) {
        print "Could not get process id for $fp_name!\n";
        push @output_data, \@blank_row;
        next;
    }
    my $crispr_id = get_crispr_id($process_id);
    if (! $crispr_id) {
        print "Could not get crispr id for $fp_name!\n";
        push @output_data, \@blank_row;
        next;
    }
    my ($crispr_seq, $wge_crispr_id) = get_crispr_info($crispr_id);
    if (! $wge_crispr_id) {
        $wge_crispr_id = '';
    } 
    push @output_row, $crispr_seq, $wge_crispr_id;
    my $crispr_strand = get_crispr_strand($crispr_id);
    if (! $crispr_strand) {
        $crispr_strand = '';
    } 
    push @output_row, $crispr_strand;
    my $forward_primer = get_forward_primer($process_id);
    if (! $forward_primer) {
        $forward_primer = '';
    }
    push @output_row, $forward_primer;
    my $reverse_primer = get_reverse_primer($process_id);
    if (! $reverse_primer) {
        $reverse_primer = '';
    }
    push @output_row, $reverse_primer;
    push @output_row, $miseq_plate_name;
    push @output_data, \@output_row;
    #last;
}

my $csv = Text::CSV->new ({ binary => 1, auto_diag => 1 });
open my $fh, ">:encoding(utf8)", "output.csv" or die "output.csv: $!";
$csv->say ($fh, $_) for @output_data;
close $fh or die "output.csv: $!";

sub get_parent_name {
    my $child_name = shift;
    my $plate = $model->schema->resultset('Plate')->find( { name => $child_name } );
    if ($plate) {
        my $parents = $plate->parent_names;
        return $$parents[0]->{name};
    }
    return;
}

sub get_miseq_exp_id {
    my $epd_name = shift;
    my $miseq_experiment = $model->schema->resultset('MiseqExperiment')->single( { name => $epd_name } );
    if ($miseq_experiment) {
        return $miseq_experiment->id;
    }
    return;
}

sub get_miseq_plate_name {
    my $miseq_exp_id = shift;
    my $miseq_id = $model->schema->resultset('MiseqExperiment')->find( { id => $miseq_exp_id } )->miseq_id;
    if ($miseq_id) {
        my $plate_id = $model->schema->resultset('MiseqPlate')->find( { id => $miseq_id } )->plate_id;
        return $model->schema->resultset('Plate')->find( { id => $plate_id } )->name;
    }
    return;
}

sub get_ep_id {
    my $ep_name = shift;
    my $ep_plate = $model->schema->resultset('Plate')->find( { name => $ep_name } );
    if ($ep_plate) {
        return $ep_plate->id;
    }
    return;
}

sub get_well_id {
    my $ep_id = shift;
    my $well_name = shift;
    my $well = $model->schema->resultset('Well')->find( { plate_id => $ep_id, name => $well_name } );
    if ($well) {
        return $well->id;
    }
    return;
}

sub get_process_id {
    my $well_id = shift;
    my $process_output_well = $model->schema->resultset('ProcessOutputWell')->find( { well_id => $well_id } );
    if ($process_output_well) {
        return $process_output_well->process_id;
    }
    return;
}

sub get_crispr_id {
    my $process_id = shift;
    my $process_crispr = $model->schema->resultset('ProcessCrispr')->find( { process_id => $process_id } );
    if ($process_crispr) {
        return $process_crispr->crispr_id;
    }
    return;
}

sub get_crispr_info {
    my $crispr_id = shift;
    my $crispr = $model->schema->resultset('Crispr')->find( { id => $crispr_id } );
    return ($crispr->seq, $crispr->wge_crispr_id);
}

sub get_crispr_strand {
    my $crispr_id = shift;
    my $crispr_locus = $model->schema->resultset('CrisprLocus')->find( { crispr_id => $crispr_id } );
    if ($crispr_locus) {
        return $crispr_locus->chr_strand;
    }
    return;
}

sub get_forward_primer {
    my $process_id = shift;
    my $design_id = $model->schema->resultset('ProcessDesign')->find( { process_id => $process_id } )->design_id;
    my $design_oligo = $model->schema->resultset('DesignOligo')->find( { design_id => $design_id, design_oligo_type_id => 'EXF' } );
    if ($design_oligo) {
        return $design_oligo->seq;
    }
    return;
}

sub get_reverse_primer {
    my $process_id = shift;
    my $design_id = $model->schema->resultset('ProcessDesign')->find( { process_id => $process_id } )->design_id;
    my $design_oligo = $model->schema->resultset('DesignOligo')->find( { design_id => $design_id, design_oligo_type_id => 'EXR' } );
    if ($design_oligo) {
        return $design_oligo->seq;
    }
    return;
}

