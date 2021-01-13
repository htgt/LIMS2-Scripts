#!/usr/bin/env perl

use strict;
use warnings;

use File::Slurp;
use Text::CSV;
use LIMS2::Model;

my $model = LIMS2::Model->new( user => 'lims2' );

my @output_data; 
my @headers = ('Gene', 'Label ID', 'Clone', 'Miseq Plate', 'Forward Primer', 'Reverse Primer', 'CRISPR Sequence', 'WGE ID', 'Strand');
push @output_data, \@headers;
my $file = $ARGV[0];
my @input_data = read_file($file, {chomp => 1});
my $counter = 0;
foreach my $line (@input_data) {
    $counter += 1;
    print "Working on line $counter\n";
    my @output_row = ();
    my ($gene, $fp_id, $clone) = split ',', $line;
    ($gene) = $gene =~ m/([A-Z0-9]+)/;
    ($fp_id) = $fp_id =~ m/([A-Za-z0-9_]+)/;
    my ($well) = $fp_id =~ m/^[A-Z0-9]+_([0-9]+)[AB_]*$/;
    $well = $well || 1;
    if ($well < 10) {
        $well = "0$well";
    }
    my $well_name = "A$well";
    ($clone) = $clone =~ m/([A-Z0-9]+)/;
    push @output_row, $gene, $fp_id, $clone;
    my $fp_name = get_fp_name($fp_id);
    if (! $fp_name) {
        print "Could not find $fp_id in database!\n";
        push @output_data, \@output_row;
        next;
    }
    my $miseq_plate_name = '';
    my $edq_name = get_edq_name($fp_name, $clone);
    if ($edq_name) {
        my $miseq_exp_id = get_miseq_exp_id($edq_name, $gene);
        if ($miseq_exp_id) {
            $miseq_plate_name = get_miseq_plate_name($miseq_exp_id);
        }
    }
    push @output_row, $miseq_plate_name;
    my $epd_name = get_parent_name($fp_name);
    if (! $epd_name) {
        print "Could not get epd name for $fp_name!\n";
        push @output_data, \@output_row;
        next;
    }
    my $ep_name = get_parent_name($epd_name);
    if (! $ep_name) {
        print "Could not get ep name for $fp_name!\n";
        push @output_data, \@output_row;
        next;
    }
    my $ep_id = get_ep_id($ep_name);
    if (! $ep_id) {
        print "Could not get ep id for $fp_name!\n";
        push @output_data, \@output_row;
        next;
    }
    my $well_id = get_well_id($ep_id, $well_name);
    if (! $well_id) {
        print "Could not get well id for $fp_name, $well_name!\n";
        push @output_data, \@output_row;
        next;
    }
    my $process_id = get_process_id($well_id);
    if (! $process_id) {
        print "Could not get process id for $fp_name!\n";
        push @output_data, \@output_row;
        next;
    }
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
    my $crispr_id = get_crispr_id($process_id);
    if (! $crispr_id) {
        print "Could not get crispr id for $fp_name!\n";
        push @output_data, \@output_row;
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
    push @output_data, \@output_row;
}

my $csv = Text::CSV->new ({ binary => 1, auto_diag => 1 });
open my $fh, ">:encoding(utf8)", "output.csv" or die "output.csv: $!";
$csv->say ($fh, $_) for @output_data;
close $fh or die "output.csv: $!";

sub get_fp_name {
    my $fp_name = shift;
	my $fp = $model->schema->resultset('Plate')->find( { name => $fp_name } );
    if (! $fp) {
        $fp_name =~ s/_/A/;
        $fp = $model->schema->resultset('Plate')->find( { name => $fp_name } );
        if (! $fp) {
            return 0;
        }
    }
    return $fp_name;
}

sub get_edq_name {
	my $fp_parent = shift;
    my $well_name = shift;
	my $plate = $model->schema->resultset('Plate')->find( { name => $fp_parent } );
	if ($plate) {
        foreach my $well ($plate->wells) {
            if ($well->well_name eq $well_name) {
                my $edq_well = $well->descendant_piq;
                if ($edq_well) {
                    return $edq_well->plate_name;
                }
            }
        }
	}
	return;
}

sub get_miseq_exp_id {
    my $edq_name = shift;
    my $gene = shift;
    my $miseq_experiment = $model->schema->resultset('MiseqExperiment')->find( { name => { like => "$edq_name%$gene%" } } );
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

sub get_parent_name {
	my $child_name = shift;
	my $plate = $model->schema->resultset('Plate')->find( { name => $child_name } );
	if ($plate) {
		my $parents = $plate->parent_names;
		return $$parents[0]->{name};
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

sub get_forward_primer {
    my $process_id = shift;
    my $process_design = $model->schema->resultset('ProcessDesign')->find( { process_id => $process_id } );
    if ($process_design) {
        my $design_oligo = $model->schema->resultset('DesignOligo')->find( { design_id => $process_design->design_id, design_oligo_type_id => 'EXF' } );
        if ($design_oligo) {
            return $design_oligo->seq;
        }
    }
    return;
}

sub get_reverse_primer {
    my $process_id = shift;
    my $process_design = $model->schema->resultset('ProcessDesign')->find( { process_id => $process_id } );
    if ($process_design) {
        my $design_oligo = $model->schema->resultset('DesignOligo')->find( { design_id => $process_design->design_id, design_oligo_type_id => 'EXR' } );
        if ($design_oligo) {
            return $design_oligo->seq;
        }
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

