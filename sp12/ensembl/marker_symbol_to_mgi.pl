#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Perl6::Slurp;
use Log::Log4perl ':easy';
use Getopt::Long;
use LIMS2::Model::Util::GeneSearch qw( retrieve_solr_gene );
use LIMS2::Util::EnsEMBL;
use Try::Tiny;
use List::MoreUtils qw( uniq );
use feature qw( say );

my $log_level = $WARN;
GetOptions(
    debug   => sub { $log_level = $DEBUG },
    verbose => sub { $log_level = $INFO },
    'species=s' => \my $species,
) and @ARGV == 1
    or die "Usage: $0 FILE\n";

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );
my $model = LIMS2::Model->new( user => 'lims2' );
my $ensembl_util = LIMS2::Util::EnsEMBL->new( species => $species );

my @genes = map { chomp; $_ } slurp( $ARGV[0] );

my %mgi_genes;
for my $gene ( @genes ) {
    my $gene_data = try{ retrieve_solr_gene( $model, { search_term => $gene } ) };

    if ( $gene_data ) {
        $mgi_genes{$gene} = $gene_data->{gene_id};
    }
    else {
        DEBUG( "Unable to find $gene in solr index" );
        my $ensembl_gene = find_ensembl_gene( $gene );
        if ( $ensembl_gene ) {
            my $external_id = external_gene_id( $ensembl_gene );
            $mgi_genes{$gene} = $external_id;
        }
        else {
            $mgi_genes{$gene} = undef;
        }
    }
}

for my $marker_symbol ( keys %mgi_genes ) {
    
    my $line = $marker_symbol . ',';
    $line .= $mgi_genes{ $marker_symbol } if $mgi_genes{ $marker_symbol };
    say $line;
}


sub find_ensembl_gene {
    my $gene = shift;

    my @genes = @{ $ensembl_util->gene_adaptor->fetch_all_by_external_name( $gene ) };
    my $gene_count = @genes;
    if ( $gene_count == 0 ) {
        WARN("Can not find EnsEMBL gene for $gene");
    }
    elsif ( @genes == 1 ) {
        return shift @genes;
    }
    else {
        WARN("Too many EnsEMBL genes found for $gene");
    }

    return;
}

sub external_gene_id {
    my  $gene = shift;

    my $type = $species eq 'Human' ? 'HGNC' : $species eq 'Mouse' ? 'MGI' : undef;
    my @dbentries = @{ $gene->get_all_DBEntries( $type ) };
    my @ids = uniq map{ $_->primary_id } @dbentries;

    if ( @ids ) {
        my $id = shift @ids;
        $id = 'HGNC:' . $id if $species eq 'Human';
        return $id;
    }

    return;
}
