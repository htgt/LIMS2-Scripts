#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Try::Tiny;
use feature qw(say);
use LIMS2::Model::Util::PrimerChecks qw(repeats_between_primer_and_target);
use WebAppCommon::Util::EnsEMBL;

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );
my $ensembl_util = {
   'Human' => WebAppCommon::Util::EnsEMBL->new( species => 'Human' ),
   'Mouse' => WebAppCommon::Util::EnsEMBL->new( species => 'Mouse' ),
};

my $primers = $model->schema->resultset("CrisprPrimer")->search({ primer_name => ['SF1','SR1','DF1','DR1'] });

my @repeats_found;
foreach my $primer ($primers->all){
	my $species = $primer->get_target->species_id;
    my $info = repeats_between_primer_and_target($primer, $ensembl_util->{$species});
    if($info){
    	push @repeats_found, $info;
    }
}

my @headings = qw(primer_id primer_name is_validated is_rejected crispr_id crispr_pair_id crispr_group_id species chromosome region_start region_end repeat_length repeat sequence);
if(@repeats_found){
	open (my $fh, ">", "repeats.csv");
	print $fh join ",", @headings;
	print $fh "\n";
    foreach my $info (@repeats_found){
    	my @out = map { $info->{$_} } @headings;
    	print $fh (join ",", @out );
    	print $fh "\n";
    }
}
