#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';

use feature qw(say);
use Data::Dumper;
use Try::Tiny;
use Text::CSV;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

sub pick_tags {
    my ($manifest, $tag, $seq, $tag_groups, $req, $read_number) = @_;
    
    my $well = $manifest->{$tag};
        
    if ($well) {
        my $read = {
            trace_indices       => $seq->trace,
            seq                 => $seq->seq,
            id                  => $seq->id,
            accession_number    => $seq->accession_number,
            qual                => $seq->qual,
            desc                => $seq->desc,
        };

        unless ($tag_groups->{$well}) {
            $tag_groups->{$well} = [];
        }
        push (@{$tag_groups->{$well}}, $read);
        say "#". $read_number . " Well: " . $well . " Tag group: " . $seq->desc;
    }

    return $tag_groups;
}

GetOptions(
    'samples=s' => \my $sample,
    'fastq=s'   => \my $fastq,
    'strand=s'  => \my $strand,
    'tag=s'     => \my $req,
);

if ($strand ne 'fwd' && $strand ne 'rev') {
    die "Incorrect strand input. Use fwd or rev";
}

my $csv = Text::CSV->new();

my $manifest;

open my $fh, '<', $sample or die "$!"; 
my $headers = $csv->getline($fh);
while (my $row = $csv->getline($fh)) {
    if ($strand eq 'fwd') {
        $manifest->{$row->[3]}->{$row->[2]} = $row->[1];
    } else {
        $manifest->{$row->[2]}->{$row->[3]} = $row->[1];
    }
}
close $fh;

my $seqio = Bio::SeqIO->new(
    -file   => $fastq,
    -format => 'fastq',
);

my $tag_reg = qw/\S:([ACTG]{8}\+[ACTG]{8})/;
my $tag_groups;
my $read_number = 1;

while (my $seq = $seqio->next_seq) {
    my ($comb_tag) = ( $seq->desc =~ $tag_reg );
    if ($comb_tag) {
        my @tags = split (/\+/, $comb_tag);

        if ($strand eq 'fwd') {
            $tag_groups = pick_tags($manifest->{$tags[1]}, $tags[0], $seq, $tag_groups, $req, $read_number);
        } else {
            $tag_groups = pick_tags($manifest->{$tags[0]}, $tags[1], $seq, $tag_groups, $req, $read_number);
        }
    }
    $read_number++;
}

say "~~~~~~~~~~~~~~~~~~~~~~~~~~";
say "Completed read collection.";
say "~~~~~~~~~~~~~~~~~~~~~~~~~~";

my @wells = keys %{$tag_groups};
my $alignment = {
    fwd => "R1",
    rev => "R2",
};

foreach my $well (@wells) {
    my $seqio_obj = Bio::SeqIO->new(
        -file => '>' . $well . '_S' . $well . '_L001_' . $alignment->{$strand} . '.fastq', 
        -format => 'fastq-sanger'
    );

    my @well_data = sort @{$tag_groups->{$well}};

    foreach my $row (@well_data) {
        my $seq_obj = Bio::Seq::Quality->new(
            -qual               => $row->{qual},
            -trace_indices      => $row->{trace_indices},
            -seq                => $row->{seq},
            -id                 => $row->{id},
            -accession_number   => $row->{accession_number},
            -desc               => $row->{desc},
        );

        $seqio_obj->write_seq($seq_obj);

        say "Writing seq to well " . $well . " - " . $row->{desc};
    }
}
say "Done!";
1;
