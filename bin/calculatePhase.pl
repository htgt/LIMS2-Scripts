#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Perl6::Slurp;
use Log::Log4perl qw( :easy );
use LIMS2::Model;
use Try::Tiny;
use Getopt::Long;

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

{
    my %log4perl = (
        level  => $INFO,
        layout => '%d %p %x %m%n'
    );

    GetOptions(
        'verbose'     => sub { $log4perl{level} = $INFO },
        'commit'      => \my $commit
    );

    Log::Log4perl->easy_init( \%log4perl );

    my @designs = split( "\n", slurp($ARGV[0]) );

    #for each design
    for my $design_id (@designs){
        #get u5_oligos data from LIMS2
        my $design;
        try {
            $design = $model->c_retrieve_design({id => $design_id});
        }
        catch {
            WARN( "$_");
        };
        next unless $design;
        #set u5_oligo_loc to u5_oligo
        my $u5_oligo = $design->oligos->find({ design_oligo_type_id => 'U5' });
        unless ($u5_oligo) {
            WARN ( "No u5 oligo can be found on design $design->id" );
            next; #design
        }
        my $locus = $u5_oligo->loci->find({assembly_id => 'GRCm38'});
        unless ($locus) {
            WARN ( "No u5 loci can be found for assembly GRCm38 on design $design->id" );
            next; #design
        }
        my %u5_oligo_loc;
        $u5_oligo_loc{start} = $locus->chr_start;
        $u5_oligo_loc{end} = $locus->chr_end;
        $u5_oligo_loc{strand} = $locus->chr_strand;
        #get transcript from LIMS2
        my $transcript = $model->ensembl_transcript_adaptor('Mouse')->fetch_by_stable_id($design->target_transcript);
        unless ($transcript) {
            WARN ( "No transcript found for $design->id" );
            next; #design
        }
        # call get_phase_from_transcript_id_and_U5_oligo
        my $phase = design_phase($transcript, \%u5_oligo_loc);
        if ($design->phase ne $phase){
            $model->txn_do(
                sub {
                INFO ( "UPDATE: updated design phase from $design->phase to $phase for design $design->id" );
                $design->update( {phase => $phase} );
                $model->txn_rollback unless $commit;
                }
            );
        }
        else {
            INFO ( "NONE: phase is already correct " . $design->phase . " for design" . $design->id )
        }
    }
}

sub design_phase {
    my ( $transcript, $u5_oligo_loc ) = @_;

    my $coding_bases = 0;
    if ( $transcript->strand == 1 ) {
        my $cs = $u5_oligo_loc->{end} + 1;
        if ( $transcript->coding_region_start > $cs or $transcript->coding_region_end < $cs ){
            return -1;
        }
        for my $e ( @{ $transcript->get_all_Exons } ){
            next unless $e->coding_region_start( $transcript );
            last if $e->seq_region_start > $cs;
            if ( $e->seq_region_end < $cs ){
                $coding_bases += $e->coding_region_end( $transcript ) - $e->coding_region_start( $transcript ) + 1;
            }
            else{
                $coding_bases += $cs - $e->coding_region_start( $transcript );
            }
        }
    }
    else{
        my $ce = $u5_oligo_loc->{start} -1;
        if ( $transcript->coding_region_start > $ce or $transcript->coding_region_end < $ce ){
            return -1;
        }
        for my $e ( @{ $transcript->get_all_Exons } ){
            next unless $e->coding_region_start( $transcript );
            last if $e->coding_region_end( $transcript ) < $ce;
            if ( $e->seq_region_start > $ce ){
                $coding_bases += $e->coding_region_end( $transcript ) - $e->coding_region_start( $transcript ) + 1;
            }
            else{
                $coding_bases += $e->coding_region_end( $transcript ) - $ce;
            }
        }
    }
    return $coding_bases %3;
}
