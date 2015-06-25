#!/usr/bin/env perl
use strict;

use LIMS2::Model;
use Perl6::Slurp;
use Try::Tiny;

my $model = LIMS2::Model->new( user => 'webapp');

my %data = map { chomp; split /,/ } slurp $ARGV[0];

$model->txn_do(
    sub {
        try{

            for my $design_id ( keys %data ) {
                my $design = $model->c_retrieve_design( { id => $design_id } );
                my $design_gene = $design->genes->first;
                $design_gene->update(
                    {
                        gene_id => $data{$design_id},
                        gene_type_id => 'MGI', 
                    }
                );
                print "modified design $design_id gene to $data{$design_id}\n";
            }

            $model->txn_rollback;
        }
        catch {
            print "failed: $_\n";
            $model->txn_rollback;
        };
    }
);
