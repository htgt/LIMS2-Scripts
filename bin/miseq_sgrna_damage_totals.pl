#! /usr/bin/perl

use strict;
use warnings;

use feature qw(say);
use Getopt::Long;
use Data::Dumper;
use File::Find;
use File::Copy qw(move);
use Text::CSV;
use LIMS2::Model;
use Try::Tiny;

GetOptions(
    'miseq=s'   => \my $miseq,
    'gene=s'    => \my $gene,
    'crispr=s'  => \my $crispr,
    'dump=s'    => \my $dump,
);

sub find_children {
    my ( $base, $reg ) = @_;
    my $fh;
    opendir ($fh, $base);
    my @files = grep {/$reg/} readdir $fh;
    closedir $fh;
    return @files;
}

sub find_file {
    my ( $base, $index, $exp, $name ) = @_;

    my $nhej_files = [];
    my $wanted = sub { _wanted( $nhej_files, $name ) };
    my $dir = $base . $exp;
    find($wanted, $dir);

    return @$nhej_files[0];
}

sub _wanted {
    return if ! -e;
    my ( $nhej_files, $file_name ) = @_;

    push(@$nhej_files, $File::Find::name) if $File::Find::name=~ /$file_name/;

    return;
}

sub read_file_lines {
    my ( $fh ) = @_;

    my @data;
    while (my $row = <$fh>) {
        chomp $row;
        push(@data, join(',', split(/\t/,$row)));
    }

    return @data;
}

sub read_columns {
    my ( $model, $csv, $fh ) = @_;

    my $overview;

    while ( my $row = $csv->getline($fh)) {
        next if $. < 2;
        my @genes;
        push @genes, $row->[1];
        $overview->{$row->[0]} = \@genes;
    }

    return $overview;
}

sub file_handling {
    my $file_name = shift;

    unless($file_name) {
        return [];
    }

    print Dumper $file_name . "\n";
    my $fh;
    open ($fh, '<:encoding(UTF-8)', $file_name) or die "$!";
    my @lines = read_file_lines($fh);
    close $fh;
    
    return \@lines;
}

sub cleavage_deletion {
    my ($sgrna, $region, $start, $length) = @_;
    my $left_bp = substr $sgrna, 0, 1;
    my $right_bp = chop $sgrna;
    if ($left_bp eq '-') {
        $sgrna = find_extent_of_deletion($sgrna, $region, $start - 1, -1);
    }

    if ($right_bp eq '-') {
        $sgrna = find_extent_of_deletion($sgrna, $region, $start + $length - 1, 1);
    }
    #Check if sgrna has bases, if not find some
    my $sgrna_len = length $sgrna;
    my $sgrna_start = index($region, lc $sgrna);
    my $hyp_count = () = $sgrna =~ /\Q-/g;
    my $hyp_perc = (100 / $sgrna_len) * $hyp_count;
    if ($hyp_perc > 30) {

        my $left_edge = $sgrna_start;
        my $right_edge = $sgrna_len;
        my $first_char = substr $sgrna, 0, 1;
        if ($first_char eq '-') {
            $left_edge = $left_edge - 5;
            $right_edge = $right_edge + 5;
        }
        if ($left_edge < 0) {
            $left_edge = 0;
        }
        if ((chop $sgrna) eq '-') {
            $right_edge += 5;
        }
        if (($right_edge + $start) > length $region) {
            $right_edge = length $region;
        }
        $sgrna = substr $region, $left_edge, $right_edge;
    }

    return $sgrna;
}

sub find_extent_of_deletion {
    my ($sgrna, $region, $pos, $direction) = @_;

    my $next_char = substr($region, $pos, 1);
    if ($next_char eq '-' || $pos < 0) {
        if ($direction == 1) {
            $sgrna = $sgrna . $next_char;
        } else {
            $sgrna = $next_char . $sgrna;
        }
        $pos = $pos + $direction; 
        $sgrna = find_extent_of_deletion($sgrna, $region, $pos, $direction);
        return $sgrna;
    }

    return $sgrna;
}


my $rna_seq = $ENV{LIMS2_RNA_SEQ};

my $base = $rna_seq . $miseq . '/';
my $experiments;

for (my $i = 1; $i < 385; $i++) {
    my $reg = "S" . $i . "_exp*" . $gene . "*";
    my @files = find_children($base, $reg);
    my @selection;
    my $percentages;
    my $classes;

    foreach my $file (@files) {
        my $allele = find_file($base, $i, $file, "Alleles_frequency_table.txt");
        my $read_count = 0;
        my $totals;

        my @lines = @{ file_handling($allele) };
        unless (@lines) {
            next;
        }
        shift @lines;
        foreach my $line (@lines) {
            my @cols = split (/,/, $line);
            my $amp = lc $cols[0];
            my $ref_seq = lc $cols[1];
            my $ref_no_hyp = $ref_seq;
            $ref_no_hyp =~ s/[-]+//g;
            my $start = index($ref_no_hyp, lc $crispr);
            my $leng = length($crispr);

            my $sgrna_ref = substr $ref_seq, $start, $leng;
            my $hyp_count = $sgrna_ref =~ tr/-//;
            $leng += $hyp_count;
            
            my $sgrna = substr $amp, $start, $leng;
            $sgrna = cleavage_deletion($sgrna, $amp, $start, $leng);
            my $primer = substr $ref_seq, $start - 15, $start;

            unless ($totals->{$sgrna}) {
                $totals->{$sgrna}->{n_deleted} = $cols[5];
                $totals->{$sgrna}->{n_inserted} = $cols[6];
                $totals->{$sgrna}->{n_mutated} = $cols[7];
                $totals->{$sgrna}->{'#Reads'} = 0;
                $totals->{$sgrna}->{'seq'} = $amp;
                $totals->{$sgrna}->{'aliseq'} = $ref_seq;
                $totals->{$sgrna}->{start} = $start;
                $totals->{$sgrna}->{leng} = $leng;
            }
            $read_count += $cols[8];
            $totals->{$sgrna}->{'#Reads'} += $cols[8];
        }
        my @reads;
        foreach my $key (keys %{ $totals }) {
            my $perc = (100 / $read_count) * $totals->{$key}->{'#Reads'};
            $totals->{$key}->{'%Reads'} = $perc;
            my $row = {
                aligned_sequence    => uc $key,
                n_deleted           => $totals->{$key}->{n_deleted},
                n_inserted          => $totals->{$key}->{n_inserted},
                n_mutated           => $totals->{$key}->{n_mutated},
                'read_count'        => $totals->{$key}->{'#Reads'},
                '%reads'            => $perc,
                seq                 => $totals->{$key}->{'seq'},
                reference_sequence  => uc $totals->{$key}->{aliseq},
                start               => $totals->{$key}->{start},
                leng                => $totals->{$key}->{leng},
            };
            push (@reads, $row);
        }

        my @sorted = sort { $b->{'read_count'} <=> $a->{'read_count'} } @reads;

        if ($dump) {
            my $csv = Text::CSV->new({binary => 1, eol => $/ })
                or die "Failed to create a CSV handle: $!";

            my $file_name = $dump . '/' . $i . "_" . $gene . "_sgrna_region_frequency.csv";
            say "Writing to $file_name\n";
            
            open my $fh, ">:encoding(utf8)", $file_name 
                or die "failed to create $file_name: $!";

            my @headings = qw( aligned_sequence n_inserted n_deleted n_mutated read_count %reads );
            $csv->column_names(\@headings);
            $csv->print($fh, [@headings]);
            foreach my $row (@sorted) {
                $csv->print ($fh, [ map { $row->{$_} } $csv->column_names ]);
            }
            
            close $fh or die "failed to close $file_name: $!";
        } else {
            print Dumper @sorted;
        }
        print "\n~~~~~~~~~~~\n";
    }
}

1;
