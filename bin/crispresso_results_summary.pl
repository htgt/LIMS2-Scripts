#! /usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;
use Path::Class;
use Getopt::Long;
use Text::CSV;

# Inputs:
#   samples: csv file containing 5 columns: Experiment ID, crispr seq, crispr strand, amplicon seq, barcode range
#   dir    : directory containing CRISPResso results directories as created by bsub_crispresso_jobs.pl
#   output : name of file to write summary info to

GetOptions(
    'samples=s' => \my $samples_file_name,
    'dir=s'     => \my $dir_name,
    'output=s'  => \my $output_name,
);

my $dir = dir($dir_name);
my $out_fh = file($output_name)->openw or die "Could not open file $output_name for writing";
my @headers = qw(experiment barcode total_reads nhej_reads nhej_percent unmodified_reads unmodified_percent output_dir);
print $out_fh (join "\t", @headers);
print $out_fh "\n";

my @samples_info = file($samples_file_name)->slurp(chomp => 1);
foreach my $line (@samples_info){
	my ($exp_id, $crispr_seq, $crispr_strand, $amplicon, $barcode_range) = split /\s*,\s*/, $line;

    my ($start,$end) = split /\s*-\s*/, $barcode_range;
    my @barcodes = $start..$end;
    foreach my $barcode (@barcodes){
    	my $bc_in_name = "_S".$barcode."_";
    	my $output_dir = "S$barcode"."_exp$exp_id";

    	my ($results_dir) = grep{ $_->basename =~ /CRISPResso_on/ } $dir->subdir($output_dir)->children;
    	if($results_dir){
            my $results_file = $results_dir->file('Alleles_frequency_table.txt');
            if($results_file){
                my $csv = Text::CSV->new({ sep_char => "\t" });
                unless(-e $results_file){
                	say "Warning: results file $results_file does not exist";
                	next;
                }
                my $fh = $results_file->openr or die "Cannot open file $results_file for reading";
                my $cols = $csv->getline($fh);
                $csv->column_names($cols);
                say "Generating summary results for $output_dir";
                my ($total, $nhej, $unmodified) = (0,0,0);
                while (my $data = $csv->getline_hr($fh)){
                	# Do the sums
                	my $num_reads = $data->{'#Reads'};
                	$total += $num_reads;
                	if($data->{NHEJ} eq 'True'){
                		$nhej += $num_reads;
                	}
                	elsif($data->{UNMODIFIED} eq 'True'){
                		$unmodified += $num_reads;
                	}
                }
                my $nhej_percent = sprintf("%.2f", ($nhej / $total) * 100);
                my $unmod_percent = sprintf("%.2f", ($nhej / $total) * 100);
                my @output = ($exp_id,$barcode,$total,$nhej,$nhej_percent,$unmodified,$unmod_percent,$output_dir);
                print $out_fh (join "\t", @output);
                print $out_fh "\n";
            }
            else{
            	say "Warning: could not find Alleles_frequency_table.txt in $results_dir";
            }
    	}
    	else{
    		say "Warning: CRISPResso results directory not found in $output_dir";
    	}
    }
}