#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Perl6::Slurp;
use Log::Log4perl ':easy';
use Getopt::Long;
use Try::Tiny;
use LIMS2::Util::EnsEMBL;
use feature qw( say );

my $log_level = $WARN;
GetOptions(
    'species=s' => \my $species,
    'gene=s' => \my $gene_id,
);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );

die( "Must set species" ) unless $species;

my $ensembl_util = LIMS2::Util::EnsEMBL->new( species => $species );

my $gene = $ensembl_util->gene_adaptor->fetch_by_stable_id( $gene_id );
unless ( $gene ) {
    WARN("Can not find EnsEMBL gene for $gene_id");
}

my $canonical_transcript = $gene->canonical_transcript;

say 'Canonical Transcript: ' . $canonical_transcript->stable_id;

my $exons = $canonical_transcript->get_all_Exons;

for my $exon ( @$exons ) {
    say 'Exon: ' . $exon->stable_id .  ' Length: ' . $exon->length;
}
