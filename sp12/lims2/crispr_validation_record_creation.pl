#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use JSON;
use Smart::Comments;

my $model = LIMS2::Model->new( user => 'webapp');

my $commit = 1;

my $crispr_es_qc_runs = $model->schema->resultset('CrisprEsQcRuns')->search(
    {
        'plate.type_id' => 'EP_PICK',
    },
    {
        join => { crispr_es_qc_wells => { well => 'plate' } },
        distinct => 1,
    }
);

while ( my $run = $crispr_es_qc_runs->next ) {
    ### qc : $run->id
    $model->txn_do(
        sub {
            update_run( $run );
            unless ( $commit ) {
                $model->txn_rollback;
            }
        }
    );
}

sub update_run {
    my $run = shift;

    for my $well ( $run->crispr_es_qc_wells->all ) {
        my $crispr = $well->crispr;
        next unless $crispr;

        my @crispr_ids;
        if ( $crispr->is_pair ) {
            @crispr_ids = ( $crispr->left_crispr_id, $crispr->right_crispr_id );
        }
        elsif ( $crispr->is_group ) {
            @crispr_ids = map { $_->crispr_id } $crispr->crispr_group_crisprs->all;
        }
        else {
            @crispr_ids = ( $crispr->id );
        }

        for my $crispr_id ( @crispr_ids ) {
            $model->schema->resultset( 'CrisprValidation' )->find_or_create(
                {
                    crispr_es_qc_well_id => $well->id,
                    crispr_id            => $crispr_id,
                }
            );
        }
    }

    return;
}
