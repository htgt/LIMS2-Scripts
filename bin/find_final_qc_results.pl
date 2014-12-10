#!/usr/bin/env perl

use strict;
use warnings;

use feature qw( say );
use Data::Dumper;

use LIMS2::Model;

my $model = LIMS2::Model->new( user => 'lims2' );

my @picks = qw(
    HFP0_Z
    HFP0001_V
    HFP0001_V_1
    HFP0001_W
    HFP0001_W_1
    HFP0001_X
    HFP0001_Y
    HFP0001_Z
    HFPS0001_V_1
    HFPS0001_V_2
    HFP0002_U_1
    HFP0002_U
    HFP0002_V
    HFP0002_X
    HFP0002_X_1
    HFP0002_Y
    HFP0003_Z
    HFP0003_Ztest
    HFP0004
    HFP0004_A_1
    HFP0004_A_3
);

say join "\t", qw(
    first_final_pick_plate
    first_final_pick_well
    parent_final_plate
    parent_final_well
    valid_final_primers
    child_final_pick_well
    design
    gene
);

for my $plate_name ( @picks ) {
    my $plate = $model->retrieve_plate( { name => $plate_name } );

    for my $well ( $plate->wells ) {

        #my $well = $plate->wells->first;

        my @parent_wells = $well->parent_wells;
        die "$plate_name is broken" if @parent_wells != 1;

        my $parent_well = $parent_wells[0];

        my $parent_plate = $parent_well->plate;

        last unless $parent_plate->type_id eq 'FINAL';

        #use gene_finder to get the marker_symbol
        my $gene = $parent_well->design->genes->first->gene_id;
        my $ms = $model->find_genes( "Human", [$gene] )->{$gene}{gene_symbol};

        my @child_wells = map { $_->output_wells } $well->child_processes;
        #die "$plate_name children are broken" if @child_wells != 1;

        my $child_fp_well = 'NO FP CHILD';
        if ( @child_wells ) {
            my $child_well = $child_wells[0];
            my $child_plate = $child_well->plate;

            if ( $plate->type_id eq 'FINAL_PICK' ) {
                $child_fp_well = $child_plate->name . "[" . $child_well->name . "]";
            }
        }
        else {
            $child_fp_well = 'NO CHILDREN';
        }

        #get genes

        # say STDERR $plate_name . " of type " . $plate->type_id;
        # say STDERR "\t" . $parent_plate->name . " of type " . $parent_plate->type;

        my $qc = $parent_well->well_qc_sequencing_result;
        #die "No qc for " . $parent_plate->name . "_" . $parent_well->name unless $qc;
        my $valid_primers = $qc ? $qc->valid_primers : '';

        say join "\t", (
            $plate_name,
            $well->name,
            $parent_plate->name,
            $parent_well->name,
            $valid_primers,
            $child_fp_well,
            $parent_well->design,
            $ms,
        );
    }
}

1;