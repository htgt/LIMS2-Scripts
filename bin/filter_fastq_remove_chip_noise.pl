#! /usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;
use Path::Class;
use feature qw/ say /;
use Text::CSV;
use Bio::Perl;
use Try::Tiny;

GetOptions(
#    'file=s'    => \my $file,
    'dir=s'     => \my $dir_name,
);
#my $csv = Text::CSV->new();

my $dir = dir($dir_name);
=head
my $params;
open my $fh, '<', $file or die "$!"; 
my $headers = $csv->getline($fh);
$csv->column_names( @{ $headers} );

while (my $row = $csv->getline_hr($fh)) {
    my $experiment = $row->{'Amplicon Name'} . '_' . $row->{'Sample Number'};
    my $fwd = $row->{'Forward primer'};
    $fwd =~ s/[^0-9]//g;
    my $well = (($fwd - 1) * 96) + $row->{'Reverse primer'};
    $params->{$well}->{$experiment} = $row->{'Amplicon'};
}
close $fh;

my @files = $dir->children;
mkdir $dir_name . '/filtered/';

my %index;
@index{@files} = (0..$#files);

foreach my $file (@files) {
    if (index($file->{file}, 'fastq') == -1) {
        my $pos = $index{$file};
        splice @files, $pos, 1;
    }
}

@files = sort { $a->{file} cmp $b->{file} } @files;

=cut

my $data;

my @files;
push @files, {file => '75_S75_L001_R1_001.fastq'};
push @files, {file => '75_S75_L001_R2_001.fastq'};
my @junk_arr = ['TTTCCCTACACGACGC','AGATCGGAAGAGCGGT','CCCGATCT','TCCGATCT','AGATCGGA','TCCCTACACGACGCTCTTCCGATCT','AGATCGGAAGAGCGGTTCAGCAGGA','CTCTTCCGATCT','AGATCGGAAGAG','TCCGATCT','AGATCGGAA']; 

my %junk = map { $_ => 1 } @junk_arr;

use Data::Dumper;
my @junk_ids;

foreach my $file (@files) {
    print Dumper $file;
    my $name = $file->{file};
    my @file_details = split /\_/, $name;
    my $well = $file_details[0];

    my $strand = 'R1';
    my $in = Bio::SeqIO->new(-format => 'fastq', -file => $dir . '/' . $file->{file});
    say $dir . '/' . $file->{file};
    my $new_file = $name;
$DB::single=1;
    my @reads;
    while (my $seq = $in->next_seq) {
        my $chip_junk = 0;
        foreach my $key(keys %junk) {
            if (index($seq->seq, $junk{$key}) != -1) {
                $chip_junk = 1;
                push(@junk_ids, $seq->id);
            }
        }
        if ($chip_junk == 0) {
            my $read = {
                trace_indices       => $seq->trace,
                seq                 => $seq->seq,
                id                  => $seq->id,
                accession_number    => $seq->accession_number,
                qual                => $seq->qual,
                desc                => $seq->desc,
            };

            #push (@{$data->{$new_file}}, $read);
            push(@reads, $read);
        }
    }
    say "$well complete.";
    my $seqio_obj = Bio::SeqIO->new(
        -file => '>' . $dir . '/filtered/' . $new_file, 
        -format => 'fastq-sanger'
    );

    say "Writing to file - $dir_name/filtered/$new_file";
    foreach my $row (@reads) {
        my $seq_obj = Bio::Seq::Quality->new(
            -qual               => $row->{qual},
            -trace_indices      => $row->{trace_indices},
            -seq                => $row->{seq},
            -id                 => $row->{id},
            -accession_number   => $row->{accession_number},
            -desc               => $row->{desc},
        );

        $seqio_obj->write_seq($seq_obj);
    }
    say "Finished writing to file - $dir_name/filtered/$well";

    
}

say "Complete";

1;
