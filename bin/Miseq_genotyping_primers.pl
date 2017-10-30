#! /usr/bin/perl

use LIMS2::Model;
use strict;
use warnings;
use Try::Tiny;
use Carp;
use Getopt::Long;
use Data::Dumper;

use WGE::Model::DB;
use LIMS2::Model::Util::OligoSelection qw(
        pick_crispr_primers
        pick_single_crispr_primers
        pick_miseq_internal_crispr_primers
        oligo_for_single_crispr
);

use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);
my $logger = Log::Log4perl->get_logger('primer_design');

my ($crispr, $design, $species);

GetOptions(
    'crispr=s'  => \$crispr,
    'species=s' => \$species,
)
or die usage_message();;

sub usage_message {
return << "END_DIE";
Usage: perl Miseq_genotyping_primers.pl
    --crispr=WGE crispr ID

END_DIE
}

my $wge_model = WGE::Model::DB->new;
my $lims2_model = LIMS2::Model->new( user => 'lims2' );
my $crispr_details = $lims2_model->schema->resultset('Crispr')->search({ id => $crispr })->first->as_hash;

$ENV{'LIMS2_SEQ_SEARCH_FIELD'} = "297";
$ENV{'LIMS2_SEQ_DEAD_FIELD'} = "0";

#my ($crispr_results2, $crispr_primers2, $chr_strand2, $chr_seq_start2) = pick_miseq_internal_crispr_primers(
my ($internal_crispr, $internal_crispr_primers) = pick_miseq_internal_crispr_primers(
    $lims2_model, {
    design_id => $design,
    crispr_id => $crispr,
    species => $species,
    repeat_mask => [''],
    offset => 20,
});

$DB::single=1;
$ENV{'LIMS2_SEQ_SEARCH_FIELD'} = "400";
$ENV{'LIMS2_SEQ_DEAD_FIELD'} = "200";

my ($pcr_crispr, $pcr_crispr_primers) = pick_miseq_internal_crispr_primers(
    $lims2_model, {
    design_id => $design,
    crispr_id => $crispr,
    species => $species,
    repeat_mask => [''],
    offset => 20,
});



$DB::single=1;

1;
