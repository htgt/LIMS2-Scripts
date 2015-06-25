#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Try::Tiny;
use feature qw( say );

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

my $plate_name = $ARGV[0];

my $plate = $model->retrieve_plate( { name => $plate_name } );

for my $well ( $plate->wells ) {
    my $crispr = try{ crispr_for_well( $well ) }; 
}

sub crispr_for_well {
    my ( $well ) = @_;

    my ( $left_crispr_well, $right_crispr_well ) = $well->left_and_right_crispr_wells;

    if ( $left_crispr_well && $right_crispr_well ) {
        my $left_crispr  = $left_crispr_well->crispr;
        my $right_crispr = $right_crispr_well->crispr;

        my $crispr_pair = $model->schema->resultset('CrisprPair')->find(
            {
                left_crispr_id  => $left_crispr->id,
                right_crispr_id => $right_crispr->id,
            }
        );

        unless ( $crispr_pair ) {
            say "$well: Unable to find crispr pair: left crispr $left_crispr, right crispr $right_crispr";
        }

        return $crispr_pair;
    }
    elsif ( $left_crispr_well ) {
        my $crispr = $left_crispr_well->crispr;
        return $crispr;
    }
    else {
        die ( "Unable to determine crispr pair or crispr for well $well" );
    }

    return;
}
