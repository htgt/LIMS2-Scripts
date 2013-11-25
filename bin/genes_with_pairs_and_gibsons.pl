#!/usr/bin/env perl

use strict;
use warnings;

use feature qw( say );

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my ( $species, @genes );
my ( $build, $compact, $pairs_only, $designs_only ) = ( 73, 0, 0, 0 );
GetOptions(
    'help'            => sub { pod2usage( 1 ) },
    'man'             => sub { pod2usage( 2 ) },
    'species=s'       => sub { my ( $name, $val ) = @_; $species = ucfirst( lc $val ) },
#    'assembly=s'      => \$assembly,
    'build=s'         => \$build,
    'genes=s{,}'      => \@genes,
    'compact!'        => \$compact,
    'pairs!'          => \$pairs_only,
    'designs!'        => \$designs_only,
) or pod2usage( 2 );

pod2usage( 1 ) unless $species and @genes;

use LIMS2::Model;
use LIMS2::Model::Util::DesignTargets qw( design_target_report_for_genes );
use List::Util qw( first );

my $model = LIMS2::Model->new( user => 'lims2' );

#we should use getopt
my %valid_species = (
    Mouse => [ 'GRCm38' ],
    Human => [ 'GRCh37' ]
);

die "Invalid species, must be one of: " . join( ", ", keys %valid_species ) 
    unless defined $valid_species{ $species };
#die "Invalid assembly, must be one of: " . join ( ", ", @{ $valid_species{$species} } )
#    unless first { $assembly eq $_ } @{ $valid_species{$species} };


my $report_params = { 
    type                 => 'simple',
    off_target_algorithm => 'bwa',
    crispr_types         => 'pair'
};

#mimic the form data (space separated gene list)
my $genes = join " ", @genes;

my ( $data ) = design_target_report_for_genes( $model->schema, $genes, $species, $build, $report_params );

# my @fields = qw( 
#     marker_symbol 
#     gene_id 
#     target_exon 
#     chromosome 
#     exon_size 
#     exon_rank 
#     designs 
#     crisprs 
#     crispr_pairs 
# );

#find all the genes with both
my ( %rows_with, %rest );
for my $row ( @{ $data } ) {
    my $gene = $row->{marker_symbol};
    if ( $row->{designs} > 0 && $row->{crispr_pairs} > 0 ) {
        push @{ $rows_with{$gene} }, $row->{ensembl_exon_id};
    }
    else {
        #allow filtering of what is displayed
        next if $pairs_only && $row->{crispr_pairs} > 0;
        next if $designs_only && $row->{designs} > 0;

        push @{ $rest{$gene} }, $row;
    }
}

#now display the genes who do
say "The following genes have both gibsons and pairs: ";
say $_ . "(" . join( ", ", @{ $rows_with{$_} } ) . ")" for keys %rows_with;

#now display all the genes that dont.
#genes will probably be in both lists because there are multiple targets per gene usually
say "The following don't have both (numbers are total_designs, total_pairs):";
while ( my ( $gene, $rows ) = each %rest ) {
    next if exists $rows_with{ $gene }; #we've already shown these

    print "$gene ";
    if ( $compact ) {
        say join " ", map { $_->{ensembl_exon_id} } @{ $rows };
    }
    else {
        say ""; #we want this on a new line
        for my $row ( @{ $rows } ) {
            say "\t" . join " ", $row->{ensembl_exon_id}, $row->{designs}, $row->{crispr_pairs};
        }
        say ""; #add a new line after each gene
    }
}

1;

__END__

=head1 NAME

genes_with_pairs_and_gibsons.pl - find all matching genes with both pairs and gibsons

=head1 SYNOPSIS

genes_with_pairs_and_gibsons.pl [options]
               
    --species          The species to fetch the data for
    --assembly         The assembly to fetch the data for
    --build            The ensembl build to search on
    --genes            One or more genes to search on
    --compact          Show just 1 line per gene
    --pairs            Only show entries that need pairs
    --designs          Only show entries that need designs
    --help             show this dialog

Example usage:

perl ./bin/genes_with_pairs_and_gibsons.pl --species human --genes CRADD BRCA2
perl ./bin/genes_with_pairs_and_gibsons.pl --species mouse --build 70 --genes Cbx1
perl ./bin/genes_with_pairs_and_gibsons.pl --species mouse --genes Cbx1 --compact --pairs
perl ./bin/genes_with_pairs_and_gibsons.pl --species Human --genes CDK5RAP2 ZNF521 LHX2 TARDBP CNOT4 --designs

=head1 DESCRIPTION

Outputs all the genes that have valid gibsons and pairs, as well as listing those that do not.

=head AUTHOR

Alex Hodgkins

=cut