#!/usr/bin/perl
use strict;
use warnings;

use LIMS2::Model;
use LIMS2::Model::Util::EngSeqParams qw( fetch_design_eng_seq_params );
use JSON qw( decode_json );
use Getopt::Long;
use Try::Tiny;
use Log::Log4perl ':easy';

my $log_level = $INFO;
GetOptions(
    'species=s'  => \my $species,
    'commit'     => \my $commit,
    'assembly=s' => \my $assembly,
);

LOGDIE( 'Must set --species' ) unless $species;
LOGDIE( 'Must set --assembly' ) unless $assembly;

my $model = LIMS2::Model->new( user => 'tasks' );
Log::Log4perl->easy_init( { level => $log_level, layout => '%p %m%n' } );

WARN ("WARNING: Only use this script after design oligos have been lifted over AND the species default assembly has been set to the new assembly");

my @qc_templates = $model->schema->resultset('QcTemplate')->search( { species_id => $species } );

MAIN: for my $qc_template ( @qc_templates ) {
    INFO( 'QC Template: ' . $qc_template->name );

    for my $qc_well ( $qc_template->qc_template_wells ) {
        my $qc_eng_seq = $qc_well->qc_eng_seq;
        my %current_params = %{ decode_json( $qc_eng_seq->params ) };
        unless ( exists $current_params{design_id} ) {
            WARN( '... QC Template plate is for crispr vectors, skipping' );
            next MAIN;
        }

        my $design = $model->c_retrieve_design( { id => $current_params{design_id} } );
        my %new_params = %{ fetch_design_eng_seq_params( $design ) };
        delete $new_params{design_type};
        delete $new_params{design_cassette_first};
        my %merged_params = %current_params;
        @merged_params{keys %new_params} = values %new_params;
        unless( $merged_params{assembly} eq $assembly ) {
            LOGDIE( "Merged eng seq params does not have $assembly coordinates: "
                    . $merged_params{assembly} );

        }
        update_eng_seq_params( $qc_eng_seq, \%merged_params );

    }
    INFO('.. Updated eng seq parameters');
}

sub update_eng_seq_params {
    my ( $qc_eng_seq, $merged_params ) = @_;

    my $json_params = JSON->new->utf8->canonical->encode($merged_params);

    $model->txn_do(
        sub {
            try{
                $qc_eng_seq->update( { params => $json_params } );
                unless ( $commit ) {
                    $model->txn_rollback;
                }
            }
            catch {
                ERROR( "Error updating qc_eng_seq record: $_" );
                $model->txn_rollback;
            };
        }
    );

}

