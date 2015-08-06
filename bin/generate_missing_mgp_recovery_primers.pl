#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Try::Tiny;
use feature qw(say);
use LIMS2::Util::QcPrimers;

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

my $primer_util = LIMS2::Util::QcPrimers->new({
    primer_project_name => 'mgp_recovery',
    model               => $model,
    base_dir            => "/nfs/users/nfs_a/af11/tmp/mgp_primers",
    persist_primers     => 1,
    overwrite           => 0,
    run_on_farm         => 0,
});

open (my $fh, ">", "/nfs/users/nfs_a/af11/tmp/mgp_primers/log.txt") or die $!;

# For each mgp repovery project
# get all experiments, get crispr groups
# if group does not have DF,DR or ER primers then generate and persis them using
# LIMS2::Util::QcPrimers
my $mgp = $model->retrieve_sponsor({ id => 'MGP Recovery'});
my @experiments = map { $_->experiments->all } $mgp->projects->all;

foreach my $exp (@experiments){
    if(my $group = $exp->crispr_group){
    	my @primer_names = map { $_->primer_name->primer_name } $group->crispr_primers;
    	if(grep { $_ eq 'DF1' } @primer_names){
            say "Crispr group ".$group->id." already has DF1 primer. Skipping.";
    	}
    	else{
    		say "Generating primers for crispr group ".$group->id;
    		### TEST THIS. used the wrong method before and did not get internal primers
            my ($picked_primers, $seq, $db_primers) = $primer_util->crispr_group_genotyping_primers($group);
            if($picked_primers){
            	print $fh $group->id.": DONE\n";
            }
            else{
            	print $fh $group->id.": FAILED\n";
            }
    	}
    }
}