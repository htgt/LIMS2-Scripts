#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Perl6::Slurp;
use Log::Log4perl ':easy';
use Getopt::Long;
use Try::Tiny;
use LIMS2::Util::EnsEMBL;
use Perl6::Slurp;
use feature qw( say );

my $log_level = $WARN;
GetOptions(
    debug   => sub { $log_level = $DEBUG },
    verbose => sub { $log_level = $INFO },
    'species=s' => \my $species,
);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %m%n' } );

die( "Must set species: --species" ) unless $species;

my $ensembl_util = LIMS2::Util::EnsEMBL->new( species => $species );

my @genes = map{ chomp; $_ } slurp \*STDIN;

my %ensembl_genes;

for my $gene ( @genes ) {
    my @genes = @{ $ensembl_util->gene_adaptor->fetch_all_by_external_name( $gene ) };
    my $gene_count = @genes;
    if ( $gene_count == 0 ) {
        WARN("Can not find EnsEMBL gene for $gene");
    }
    elsif ( @genes == 1 ) {
        $ensembl_genes{$gene} = $genes[0]->stable_id;
    }
    else {
        WARN("Too many EnsEMBL genes found for $gene");
    }

}

for my $gene ( @genes ) {
    my $ens_gene = exists $ensembl_genes{ $gene } ? $ensembl_genes{ $gene } : '';
    my $line = $gene . ',' . $ens_gene; 
    say $line;
}
