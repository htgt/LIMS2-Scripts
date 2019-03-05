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

sub find_children {
    my ( $base, $reg ) = @_;
    my $fh;
    opendir ($fh, $base);
    my @files = grep {/$reg/} readdir $fh;
    closedir $fh;
    return @files;
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

sub file_handling {
    my $file_name = shift;

    my $fh;
    open ($fh, '<:encoding(UTF-8)', $file_name) or die "$!";
    my @lines = read_file_lines($fh);
    close $fh;
    
    return \@lines;
}

sub find_file {
    my ( $base, $index, $exp, $name ) = @_;

    my $nhej_files = [];
    my $wanted = sub { _wanted( $nhej_files, $name ) };
    my $dir = $base . "S" . $index . "_exp" . $exp;
    find($wanted, $dir);

    return @$nhej_files[0];
}

sub _wanted {
    return if ! -e;
    my ( $nhej_files, $file_name ) = @_;

    push(@$nhej_files, $File::Find::name) if $File::Find::name=~ /$file_name/;

    return;
}

my $rna_seq = $ENV{LIMS2_RNA_SEQ};

my $folder_reg = "Miseq_[0-9]+.*";
my @folders = find_children($rna_seq, $folder_reg);

my $experiments;
foreach my $folder (@folders) {
    my $base = $rna_seq . $folder . '/';
    print "Analysing " . $folder . "\n";
    for (my $i = 1; $i < 385; $i++) {
        my $reg = "S" . $i . "_exp[A-Za-z0-9_]+";
        my @files = find_children($base, $reg);
        my @exps;
        foreach my $file (@files) {
            my @matches = ($file =~ /S\d+_exp([A-Za-z0-9_]+)/g);
            foreach my $match (@matches) {
                push (@exps,$match);
            }
        }

        @exps = sort @exps;
        my @selection;

        foreach my $exp (@exps) {
            my $quant = find_file($base, $i, $exp, "Quantification_of_editing_frequency.txt");

            if ($quant) {
                my @lines = @{file_handling($quant)};
                my @wt = ($lines[1] =~ qr/^,- Unmodified:(\d+)/);
                my @nhej = ($lines[2] =~ qr/^,- NHEJ:(\d+)/);
                my @hdr = ($lines[3] =~ qr/^,- HDR:(\d+)/);
                my @total = ($lines[6] =~ qr/^Total Aligned:(\d+)/);

                $experiments->{$folder}->{$exp}->{sprintf("%02d", $i)} = {
                    wt      => $wt[0],
                    nhej    => $nhej[0],
                    hdr     => $hdr[0],
                    total   => $total[0],
                };
            }
        }
    }
}
my $result;
my $run_ids;
my $gene_exp;

my $model = LIMS2::Model->new({ user => 'tasks' });
foreach my $folder (keys %{$experiments}) {
    my $miseq_db = $model->schema->resultset('Plate')->find({ name => $folder })->miseq_details;
    $run_ids->{$folder} = $model->schema->resultset('MiseqPlate')->find({ id => $miseq_db->{id} })->as_hash->{run_id};
    foreach my $exp (keys %{$experiments->{$folder}}) {
        my $nhej = 0;
        my $hdr = 0;
        my $wt = 0;
        my $total = 0;
        foreach my $index (keys %{$experiments->{$folder}->{$exp}}) {
            $nhej += $experiments->{$folder}->{$exp}->{$index}->{nhej};
            $hdr += $experiments->{$folder}->{$exp}->{$index}->{hdr};
            $wt += $experiments->{$folder}->{$exp}->{$index}->{wt};
            $total += $experiments->{$folder}->{$exp}->{$index}->{total};
        }
        if ($total == 0) {
            $total++;
        }
        my $nhej_target = sprintf("%.2f", ($nhej / $total) * 100);
        my $hdr_target = sprintf("%.2f", ($hdr / $total) * 100);
        my $wt_target = sprintf("%.2f", ($wt / $total) * 100);
        $result->{$folder}->{$exp} = {
            nhej    => $nhej,
            hdr     => $hdr,
            wt      => $wt,
            total   => $total,
            nhej_efficiency => $nhej_target,
            hdr_efficiency => $hdr_target,
            wt_efficiency => $wt_target,
        };
        print "Plate: " . $folder . ", Experiment: " . $exp .", WT: " . $wt . ", NHEJ: " . $nhej. ", HDR: " . $hdr . ", Total: " . $total . ", NHEJ_Eff: " . $nhej_target . "%, HDR_Eff: " . $hdr_target . "%\n";
        my $comp_key = $folder . '_' . $exp;
        try {
            $gene_exp->{$comp_key} = $model->schema->resultset('MiseqExperiment')->find({ miseq_id => $miseq_db->{id}, name => $exp })->gene;
        };
    }

}

use Data::Dumper;
print Dumper $run_ids;
print Dumper $gene_exp;

my @headers = qw( Run_ID Miseq Experiment Gene WT_reads NHEJ_reads HDR_reads Total_reads WT% NHEJ% HDR% );
my $csv = Text::CSV->new({binary => 1, eol => "\n"}) or die "Cannot use CSV: ".Text::CSV->error_diag ();
$csv->column_names(\@headers);
my $new_file = 'NHEJ_HDR_targeting_efficiencies.csv';
open my $fh, '>', $new_file or die "$new_file: $!";
$csv->print ($fh, \@headers);
my @keys = sort keys %{$result};
foreach my $folder_name (@keys) {
    print "Printing $folder_name \n";
    my $miseq = $result->{$folder_name};
    my $run_id = $run_ids->{$folder_name};
    foreach my $experiment (keys %{$miseq}) {
        my $exp_details = $miseq->{$experiment};
        my $row = {
            Run_ID          => $run_id,
            Miseq           => $folder_name,
            Experiment      => $experiment,
            Gene            => $gene_exp->{$folder_name . '_' . $experiment},
            WT_reads        => $exp_details->{wt},
            NHEJ_reads      => $exp_details->{nhej},
            HDR_reads       => $exp_details->{hdr},
            Total_reads     => $exp_details->{total},
            'WT%'           => $exp_details->{wt_efficiency} . "%",
            'NHEJ%'         => $exp_details->{nhej_efficiency} . "%",
            'HDR%'          => $exp_details->{hdr_efficiency} . "%",
        };
        #$csv->print_hr($fh, $row);
        $csv->print ($fh, [ map { $row->{$_} } $csv->column_names ]);
    }
}
close $fh;
print "Finished.";
