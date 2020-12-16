#!/usr/bin/env perl

use strict;
use warnings;

use Text::CSV;
use LIMS2::Model;

=pod

This script was written to retrieve total reads and NHEJ reads for each well in the arrayed plates on Miseq 96

=cut

my $model = LIMS2::Model->new( user => 'lims2' );

my @headers = ( 'Well name', 'Total reads', 'NHEJ reads' );

foreach my $experiment ( 'ES4_1', 'ES4_2', 'ES4_3', 'HCC_1', 'HCC_2', 'HCC_3' )
{
    my @output_data      = ( \@headers );
    my $miseq_experiment = $model->schema->resultset('MiseqExperiment')
      ->find( { name => $experiment } );
    if ( !defined $miseq_experiment ) {
        warn "$experiment not found in database";
        next;
    }
    my @miseq_well_experiments =
      map { $_->as_hash }
      $model->schema->resultset('MiseqWellExperiment')
      ->search( { miseq_exp_id => $miseq_experiment->id } );

    foreach my $miseq_well_exp (@miseq_well_experiments) {
        my $well_name = $model->schema->resultset('Well')
          ->find( { id => $miseq_well_exp->{well_id} } )->well_name;
        my @output_row = (
            $well_name,
            $miseq_well_exp->{total_reads},
            $miseq_well_exp->{nhej_reads}
        );
        push @output_data, \@output_row;
    }

    my $csv = Text::CSV->new( { binary => 1, auto_diag => 1 } );
    open my $fh, ">:encoding(utf8)", "$experiment.csv"
      or die "$experiment.csv: $!";
    $csv->say( $fh, $_ ) for @output_data;
    close $fh or die "$experiment.csv: $!";
}

