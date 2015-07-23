#! /usr/bin/perl

use LIMS2::Model;
use strict;
use warnings;

use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);
my $logger = Log::Log4perl->get_logger('count_crispr_damage');

my $lims2_model = LIMS2::Model->new( { user => 'lims2' } );
$logger->info('Connecting to LIMS2');
my $lims2_schema = $lims2_model->schema;

$logger->info('Retrieving crispr validation data');
my @crispr_validations = $lims2_schema->resultset( 'CrisprValidation' )->search(
     { },
     {
        join => [ qw/ crispr crispr_es_qc_well / ],
    }
);

$logger->info( scalar(@crispr_validations) . ' rows were returned');

my %clone_clip;

foreach my $crispr_val ( @crispr_validations ) {
my %crispr_clip;
    # deal at the crispr level first
            
    if ( $crispr_val->validated == 0 ) {
        $crispr_clip{ $crispr_val->crispr_id }->{ 'invalid' } += 1;
    }
    elsif ( $crispr_val->validated == 1) {
        $crispr_clip{ $crispr_val->crispr_id }->{ 'valid' }  += 1;
    }
    if ( !exists $crispr_clip{ $crispr_val->crispr_id }->{ 'crispr_seq' }) {
        $crispr_clip{ $crispr_val->crispr_id }->{ 'crispr_seq' } = $crispr_val->crispr->seq;
    }
    #add it to the clone level clip
    $clone_clip{ $crispr_val->crispr_es_qc_well->well_id } = \%crispr_clip ;
}

$logger->info( 'Run Completed');

