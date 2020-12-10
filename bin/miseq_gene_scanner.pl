#!/usr/bin/env perl

use strict;
use warnings;

use LIMS2::Model;

my $model = LIMS2::Model->new( user => 'lims2' );

my @miseq_exps = map { $_->as_hash } $model->schema->resultset('MiseqExperiment')->search;
printf( "%25s %25s %25s\n", "Experiment", "Gene", "Solr Gene" );
foreach my $miseq_exp (@miseq_exps) {
    my $exp_name  = $miseq_exp->{name};
    my $gene_name = $miseq_exp->{gene};
    if ( $miseq_exp->{experiment_id} ) {
        my $experiment
            = $model->schema->resultset('Experiment')->find( { id => $miseq_exp->{experiment_id} } )->as_hash;
        my $solr_gene = $model->find_gene(
            {   search_term => $experiment->{gene_id},
                species     => 'Human',
            }
        );
        my $solr_gene_symbol = $solr_gene->{gene_symbol};
        if ( ( not $exp_name =~ /$solr_gene_symbol/ ) || ( not $exp_name =~ /$gene_name/ ) ) {
            printf( "%25s %25s %25s\n", $exp_name, $gene_name, $solr_gene_symbol );
        }
    }
    else {
        if ( not $exp_name =~ /$gene_name/ ) {
            printf( "%25s %25s %25s\n", $exp_name, $gene_name, "N/A" );
        }
    }
}
