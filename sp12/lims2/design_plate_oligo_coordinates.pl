#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use feature qw( say );

use Smart::Comments;

# Print design plate oligo coordinates

my $model = LIMS2::Model->new( user => 'webapp' );

my $plate = $model->retrieve_plate( { name => $ARGV[0] } );

say 'design_id,gene,well_name,chromosome,oligo_type,start,end,strand';
for my $well ( $plate->wells->all ) {
    my $design = $well->design;
    my @oligos = $design->oligos;
    my %oligo_data = map { my $d = $_->as_hash; $d->{type} => $d->{locus} } @oligos;

    my $gene = $design->genes->first->gene_id;
    my $common = $design->id . ',' . $gene . ',' . $well->name . ',' . $oligo_data{U5}{chr_name} . ',';
    for my $d ( keys %oligo_data ) {
        say $common . $d . ',' . $oligo_data{$d}{chr_start} . ',' . $oligo_data{$d}{chr_end} . ',' . $oligo_data{$d}{chr_strand};
    }
}

