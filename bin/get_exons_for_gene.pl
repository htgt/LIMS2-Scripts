#!/usr/bin/perl

use strict;
use warnings;

use LIMS2::Util::EnsEMBL;
use List::MoreUtils qw( uniq );

#given a file of params for paired_crisprs.sh
#which is in the format 'Marker_Symbol Species Exon_IDs'
#fill in any missing exons from the canonical transcript

die "Usage: get_exons_for_gene.pl species file.txt" unless @ARGV > 1

my $species = shift;

my $e = LIMS2::Util::EnsEMBL->new( species => $species );

while ( my $line = <> ) {
    my ( $marker_symbol, $row_species, @exon_ids ) = split " ", $line;

    die "Species don't match!" unless lc($species) eq lc($row_species);

    my $gene = $e->gene_adaptor->fetch_by_exon_stable_id( $exon_ids[0] );

    die "Couldn't find a gene for " . $exon_ids[0] unless $gene;

    push @exon_ids, @{$gene->canonical_transcript->get_all_Exons};

    say join ",", $marker_symbol, $row_species, uniq @exon_ids;
}