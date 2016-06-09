#!/usr/bin/env perl

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use feature qw(say);

my $input = Bio::SeqIO->new(-file   => $ARGV[0], -format => 'genbank');

my $out_name = $ARGV[0];
$out_name =~ s/gbk/reverse.gbk/;

open (my $out_fh, ">", $out_name) or die $!;

my $output = Bio::SeqIO->new(-fh => $out_fh, -format => 'genbank');

my $seq = $input->next_seq;

my @feat=$seq->get_SeqFeatures();
my $seq_length = $seq->length();
my $rev_seq = $seq->revcom;

for my $feat(@feat){
	say "-----------------------";
	say "start  : ".$feat->start;
	say "end    : ".$feat->end;
	say "length : ".$feat->length;

	my $feat_length = $feat->length();

    my $end=($seq_length - $feat->location->start) + 1;
    my $start=($end - $feat_length) + 1;

    my $strand=$feat->strand;

    $strand=$strand*-1;


    next if(ref($feat->location) ne "Bio::Location::Simple");

    $feat->location->start($start);
    $feat->location->end($end);
    $feat->strand($strand);

	say "*start  : ".$feat->start;
	say "*end    : ".$feat->end;
	say "*length : ".$feat->length;    
    $rev_seq->add_SeqFeature($feat);
}

$output->write_seq($rev_seq);
say "Created file $out_name";

