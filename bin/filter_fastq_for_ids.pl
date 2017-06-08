#! /usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;
use Path::Class;
use feature qw/ say /;

GetOptions(
    'dir=s'      => \my $dir_name,
);


$DB::single=1;
my $dir = dir($dir_name);
my @files = $dir->children;
mkdir $dir_name . '/ids/';

@files = sort { $a->{file} cmp $b->{file} } @files;

foreach my $file (@files) {
    say $file->{file};
    my $in = Bio::SeqIO->new(-format => 'fastq', -file => './' . $file->{file});
    my @ids;
    while (my $seq = $in->next_seq) {
        if (index($seq->seq, "GAATGG") != -1) {
            push @ids, $seq->display_id;
        }
    }
    my @names = split(/_/, $file->{file});
    my $file_name = $names[0];
    open (my $fh, '>', './ids/' . $file_name . '.txt');
    print $fh join ("\n", @ids);
    close $fh;
}
