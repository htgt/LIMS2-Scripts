#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use feature qw( say );

use LIMS2::Model;
use LIMS2::Model::Util::DesignInfo;

use Smart::Comments;

my $model = LIMS2::Model->new( user => 'webapp');

my $design_id = $ARGV[0];

my $design = $model->retrieve_design( { id => $design_id } );

# got design

my $di = LIMS2::Model::Util::DesignInfo->new( { design => $design } );

say 'Type: ' . $di->type;
say 'Assembly: ' . $di->default_assembly;
say 'Chromosome: ' . $di->chr_name;
say 'Strand: ' . $di->chr_strand;

say 'Homology Arm Start: ' . $di->homology_arm_start;
say 'Homology Arm End: ' . $di->homology_arm_end;

say 'Target Region Start: ' . $di->target_region_start;
say 'Target Region End: ' . $di->target_region_end;

say 'Cassette Start: ' . $di->cassette_start;
say 'Cassette End: ' . $di->cassette_end;

say 'Loxp Start: ' . $di->loxp_start;
say 'Loxp End: ' . $di->loxp_end;


say $di->target_gene->stable_id;
$di->target_gene->external_db('HGNC');
say $di->target_gene->external_name;

my $exons = $di->floxed_exons;

for my $exon ( @{ $exons } ) {
    say $exon->stable_id;
}

