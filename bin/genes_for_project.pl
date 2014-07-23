#!/usr/bin/env perl

use strict;
use warnings;

use feature qw( say );
use Data::Dumper;

use LIMS2::Model;
use LIMS2::Util::EnsEMBL;

my $model = LIMS2::Model->new( user => 'lims2' );
my $ens = LIMS2::Util::EnsEMBL->new( species => "Human" );

my @genes = $model->schema->resultset("Project")->search(
    { sponsor_id => "Transfacs", targeting_type => "single_targeted", species_id => "Human" }
);

my @cols = qw(
    marker_symbol
    ensembl_gene_id
    ensembl_exon_id
    exon_chr_name
    exon_chr_start
    exon_chr_end
);

say join ",", @cols;

for my $project_gene ( @genes ) {
    #has gene_id, design_count, gene_symbol, crispr_pairs_count, ensembl_id
    my $gene_data = $model->find_gene(
        {
            search_term => $project_gene->gene_id,
            species => "Human",
            show_all => 1
        }
    );

    my $gene = $ens->gene_adaptor->fetch_by_stable_id( $gene_data->{ ensembl_id } );

    unless ( $gene ) {
        say STDERR "Skipping " . $gene_data->{ ensembl_id } . " can't be found!" ;
        next;
    }

    for my $exon ( @{ $gene->canonical_transcript->get_all_Exons } ) {
        say join ",", $gene->external_name,
                      $gene->stable_id,
                      $exon->stable_id,
                      $exon->seq_region_name,
                      $exon->seq_region_start,
                      $exon->seq_region_end;
    }
}

__END__

=head1 NAME

genes_for_project.pl - retrieve all the different gene ids for a project as csv

=head1 SYNOPSIS

genes_for_project.pl [options]

    --help            Show this dialog

Example usage:

./genes_for_project.pl

=head1 DESCRIPTION

Create a CSV of exon data for all genes in a project

=head AUTHOR

Alex Hodgkins

=cut
