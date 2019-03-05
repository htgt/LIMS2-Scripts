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
use Data::Dumper;


sub alphanumeric_sort {
    my ($a, $b) = @_;
    
    my ($number_a) = $a =~ /(\d+)_/;
    my ($number_b) = $b =~ /(\d+)_/;

    my ($letter_a) = $a =~ /\S+R(1-2)/;
    my ($letter_b) = $b =~ /\S+R(1-2)/;

    return $number_a <=> $number_b or $letter_a <=> $letter_b;
}


GetOptions(
    'file=s'    => \my $file_name,
    'dir=s'     => \my $dir_name,
    'seq=s'     => \my $seq_segment,
    'gene=s'    => \my $gene_name,
);
my $params;
if ($file_name) {
    my $csv = Text::CSV->new();

    open my $fh, '<', $file_name or die "$!"; 
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
} elsif ($seq_segment) {
    my @segments = split /,/, $seq_segment;
    for (my $i = 1; $i < 385; $i++) { 
        $params->{$i}->{$gene_name} = [@segments];
    }
} else {
    say "No file or sequence specified";
    return;
}
my $dir = dir($dir_name);
my @files = sort { alphanumeric_sort($a->{file}, $b->{file}) } grep {/\.fastq/} $dir->children;
mkdir $dir_name . '/filtered/';
my $data;

foreach my $file (@files) {
    say "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
    my $name = $file->{file};
    say $name;
    my @file_details = split /\_/, $name;
    my $well = $file_details[0];
    unless ($params->{$well}) {
        next;
    }

    my $strand = $file_details[3];
    my $in = Bio::SeqIO->new(-format => 'fastq', -file => $dir . '/' . $file->{file});
    my $file_reads;
    foreach my $gene (keys %{$params->{$well}}) {
        say "Started $gene - $well";
        my @amps = @{ $params->{$well}->{$gene} };
        my $new_file = $gene . '_' . $name;
        my @reads;
        my $read_counts;
        while (my $seq = $in->next_seq) {
            foreach my $amp (@amps) {
                my $fwd_amp = substr ($amp, 0, 25);
                my $rev_amp = revcom($fwd_amp);
                my $amp_str = $fwd_amp;
                if ($strand eq 'R2') {
                    $amp_str = $rev_amp->seq;
                }

                if (index($seq->seq, $amp_str) != -1) {
                    unless ($read_counts->{$amp_str}) {
                        $read_counts->{$amp_str} = 0;
                    }
                    $read_counts->{$amp_str} += 1;
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
        }
        my $count = scalar @reads;
        $file_reads->{$strand} = $count;
        print Dumper $read_counts;
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
        say "\n";
        say "$gene - $well $strand complete. $count reads found.";
    }
}

say "Complete";

1;
