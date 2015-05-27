#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Log::Log4perl ':easy';
use Getopt::Long;
use Try::Tiny;
use LIMS2::Util::EnsEMBL;
use List::MoreUtils qw( uniq );

use Smart::Comments;

GetOptions(
    'species=s'            => \my $species,
    'canonical_transcript' => \my $canonical_transcript,
    'all_exons'            => \my $all_exons,
    'canonical_exons'      => \my $canonical_exons,
    'external_id'          => \my $external_id,
);

Log::Log4perl->easy_init( { level => $WARN, layout => '%p %m%n' } );

LOGDIE( 'Must set species' ) unless $species;

my $ensembl_util = LIMS2::Util::EnsEMBL->new( species => $species );

while ( my $gene_id = <> ) {
    chomp($gene_id);
    my $gene = $ensembl_util->gene_adaptor->fetch_by_stable_id( $gene_id );
    unless ( $gene ) {
        WARN("Can not find EnsEMBL gene for $gene_id");
        next;
    }
    print "$gene_id";

    my $ct = $gene->canonical_transcript;
    if ( $canonical_transcript ) {
        print  ',' . $ct->stable_id;
    }

    if ( $canonical_exons ) {
        my @exons = @{ $ct->get_all_Exons };
        print ',' .  join( '|', map { $_->stable_id } @exons );
    }

    if ( $all_exons ) {
        my @exons = @{ $gene->get_all_Exons };
        print ',' .  join( '|', map { $_->stable_id } @exons );
    }

    if ( $external_id ) {
        my $id = external_gene_id( $gene );
        print ',' . $id if $id;
    }

    print "\n";
}

sub external_gene_id {
    my ( $gene ) = @_;

    my $type = $species eq 'Human' ? 'HGNC' : $species eq 'Mouse' ? 'MGI' : undef;
    my @dbentries = @{ $gene->get_all_DBEntries( $type ) };
    my @ids = uniq map{ $_->primary_id } @dbentries;
    ### @ids

    if ( @ids ) {
        my $id = shift @ids;
        $id = 'HGNC:' . $id if $species eq 'Human';
        return $id;
    }

    return;
}
