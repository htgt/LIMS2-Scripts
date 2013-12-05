#!/usr/bin/env perl

use strict;
use warnings;

use feature qw( say );

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Try::Tiny;

my ( @exon_ids );
GetOptions(
    'help'            => sub { pod2usage( 1 ) },
    'man'             => sub { pod2usage( 2 ) },
    'exon-ids=s{,}'   => \@exon_ids,
) or pod2usage( 2 );

pod2usage( 1 ) unless @exon_ids;

use LIMS2::Model;
use List::Util qw( first );

my $model = LIMS2::Model->new( user => 'lims2' );

#get the ids to give to DesignTargetCrisprs
my @design_targets = map { $_->id } $model->schema->resultset("DesignTarget")->search(
    { ensembl_exon_id => { -IN => [ @exon_ids ] } }
);

die "No targets found" unless @design_targets;

my @data = $model->schema->resultset("DesignTargetCrisprs")->search(
    { 
        design_target_id => { -IN => [ @design_targets ] }, 
        "off_target_summaries.algorithm" => "bwa" 
    }, 
    { prefetch => { crispr => "off_target_summaries" } }
);

for my $row ( @data ) {
    my $crispr = $row->crispr;

    for my $check ( qw( crispr_designs process_crisprs ) ) {
        if ( $crispr->$check->all ) {
            die "Can't delete crispr " . $crispr->id . " as it's linked to $check";
        }
    }

    my $pairs_rs;
    if ( $crispr->pam_right ) {
        $pairs_rs = $crispr->crispr_pairs_right_crisprs;
    }
    else {
        $pairs_rs = $crispr->crispr_pairs_left_crisprs;
    }

    if ( map { $_->crispr_designs } $pairs_rs->all ) {
        die "Can't delete " . $crispr->id . " as a pair is linked to a design.";
    }

    #everything looks ok to delete so lets do it
    $model->txn_do(
        sub {
            try{
                $pairs_rs->delete;
                $crispr->off_target_summaries->delete;
                $crispr->off_targets->delete;
                $crispr->loci->delete;
                $crispr->delete;
            }
            catch {
                say "Deletion failed: $_";
                $model->txn_rollback;
            };
        }
    );
}

1;

__END__

=head1 NAME

remove_crisprs_for_exons.pl - delete all crisprs and pairs linked to an exon

=head1 SYNOPSIS

remove_crisprs_for_exons.pl [options]
               
    --exon-ids         The exons whose crisprs you want to delete
    --help             show this dialog

Example usage:

perl remove_crisprs_for_exons.pl --exon-ids ENSE00001353200 ENSE00001140991

=head1 DESCRIPTION

Given an exon id this script will delete every crispr associated (+/- 200bp), and
the related tables.

=head AUTHOR

Alex Hodgkins

=cut