#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';

use feature qw(say);
use Data::Dumper;
use TryCatch;
use LIMS2::Model::Util::RedmineAPI;

use Log::Log4perl qw( :easy );

BEGIN { Log::Log4perl->easy_init( { level => $DEBUG } ) }

# Takes a csv file containing lines of the format:
# <redmine issue ID>,<new status>

my $redmine = LIMS2::Model::Util::RedmineAPI->new_with_config();

my $file = $ARGV[0];

open (my $fh, "<", $file) or die "Could not open file $file for reading - $!";

foreach my $line (<$fh>){
	chomp $line;
	my ($id,$status) = split /\s*,\s*/, $line;
	try{
    	$redmine->update_issue_status($id,$status);
    	say "Updated issue $id to status $status";
    }
    catch($err){
    	say "Could not update issue $id with status $status: $err";
    }
}