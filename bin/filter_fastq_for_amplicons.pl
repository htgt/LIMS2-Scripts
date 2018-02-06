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
    'file=s'    => \my $file,
    'dir=s'     => \my $dir_name,
);
my $csv = Text::CSV->new();

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


my $dir = dir($dir_name);
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

my $data;

foreach my $file (@files) {
    my $name = $file->{file};
    say $name;
    my @file_details = split /\_/, $name;
    my $well = $file_details[0];
    unless ($params->{$well}) {
        next;
    }

    my $strand = $file_details[3];
    my $in = Bio::SeqIO->new(-format => 'fastq', -file => $dir . '/' . $file->{file});
    foreach my $gene (keys %{$params->{$well}}) {
        say "Started $gene - $well";
        my $amp = $params->{$well}->{$gene};
        my $fwd_amp = substr ($amp, 0, 25);
        my $rev_amp = revcom($fwd_amp);
        my $amp_str = $fwd_amp;
        if ($strand eq 'R2') {
            $amp_str = $rev_amp;
        }

        my $new_file = $gene . '_' . $name;

        my @reads;
        while (my $seq = $in->next_seq) {
            if (index($seq->seq, $amp_str) != -1) {
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
        say "$gene - $well complete.";
        my $seqio_obj = Bio::SeqIO->new(
            -file => '>' . $dir_name . '/filtered/' . $new_file, 
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
}

say "Complete";

1;
