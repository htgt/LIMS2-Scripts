#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Try::Tiny;
use WebAppCommon::Util::EnsEMBL;
use List::MoreUtils qw( uniq );
use Text::CSV;

my $ensembl_util = WebAppCommon::Util::EnsEMBL->new( species => 'Mouse' );

# ensembl_gene_id,marker_symbol
my $file = $ARGV[0];
open ( my $fh, '<', $file ) or die( "Can not open $file " . $! );
my $csv = Text::CSV->new();
$csv->column_names( @{ $csv->getline( $fh ) } );

my @COLUMN_HEADERS = qw(
ensembl_gene_id
marker_symbol
description
biotype
);

my $output = IO::Handle->new_from_fd( \*STDOUT, 'w' );
my $output_csv = Text::CSV->new( { eol => "\n" } );
$output_csv->print( $output, \@COLUMN_HEADERS );

while ( my $data = $csv->getline_hr( $fh ) ) {
    my %data;
    my $gene = $ensembl_util->gene_adaptor->fetch_by_stable_id( $data->{ensembl_gene_id} );

    $data{ensembl_gene_id} = $data->{ensembl_gene_id};
    $data{marker_symbol} = $data->{marker_symbol};
    $data{description} = $gene->description;
    $data{biotype} = $gene->biotype;

    $output_csv->print( $output, [ @data{ @COLUMN_HEADERS } ] );
}
