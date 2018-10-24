#! /usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;
use Path::Class;
use Getopt::Long;
use Text::CSV;

# Inputs:
#   samples: csv file containing 6 columns: Experiment ID, crispr seq, crispr strand, amplicon seq, minimum index and maximum index
#   dir    : directory containing input fastq files (1 per barcode with names containing _S<barcode number>_ and R1 or R2)

# NB: Experiment ID is just a unique identifier for the combination of amplicon and crispr used.
# It is not a LIMS2 experiment ID at this stage in development.

my $crispresso = $ENV{CRISPRESSO_CMD} || "/nfs/team87/farm3_lims2_vms/software/python_local/bin/CRISPResso";

GetOptions(
    'samples=s' => \my $samples_file_name,
    'dir=s'     => \my $dir_name,
    'offset=s'  => \my $offset,
    'single'    => \my $single,
);

# Ouputs:
# 1 directory per barcode/experiment combination
# directory names in format S<barcode number> _exp<experiment ID>

# bsub -q normal -G team87-grp -eo S1_expA_farm/job.err -oo S1_expA_farm/job.out -M500 -R"select[mem>500] rusage[mem=500]" ./bin/CRISPResso -w 30 --hide_mutations_outside_window_NHEJ --ignore_substitutions -o S1_expA -r1 real_data/orig_fastq/Homo-sapiens_S1_L001_R1_001.fastq -a GAAAGTCCGGAAAGACAAGGAAGgaacacctccacttacaaaagaagataagacagttgtcagacaaagccctcgaaggattaagccagttaggattattccttcttcaaaaaggacagatgcaaccattgctaagcaactcttacagag -g CTCGAAGGATTAAGCCAGTT

$offset = $offset || 0;

my $csv = Text::CSV->new({ binary => 1 }) or die "Can't use CSV: " . Text::CSV->error_diag();
open my $fh, '<:encoding(UTF-8)', $samples_file_name or die "Can't open CSV: $!";

my $dir = dir($dir_name);
my @files = $dir->children;
my @barcodes;

my @col = $csv->column_names($csv->getline($fh));
if ($col[5] =~  m/^min/gmi){
    my @heads = qw(experiment gene crispr strand amplicon min_index max_index hdr);
    $csv->column_names(\@heads);
} else{
    my @heads = qw(experiment gene crispr strand amplicon range hdr);
    $csv->column_names(\@heads);
}

while (my $hr = $csv->getline_hr($fh)){
    if (exists $hr->{range}){
        my @barcode_sets = split /\s*;\s*/, $hr->{range};
        foreach my $wells (@barcode_sets) {
            if ($wells =~ /\s*-\s*/){
                my ($start,$end) = split /\s*-\s*/, $wells;
                push (@barcodes, $start..$end);
            } else {
                push (@barcodes, $wells);
            }
        }
    }
    elsif (exists $hr->{min_index} && $hr->{max_index}){
        push (@barcodes, $hr->{min_index}..$hr->{max_index});
    }
    else{
        warn "Range is not provided in any of the expected formats.";
        die;
    }
    my $crispr_site = substr($hr->{crispr},0,20);
	foreach my $barcode (@barcodes) {

	    my $bc_in_name = "_S".$barcode."_";
        my ($fwd_file) = grep { $_=~ /$bc_in_name.*R1/ } @files;
        my ($rev_file) = grep { $_=~ /$bc_in_name.*R2/ } @files;

        my $fix = $barcode - $offset;
        my $output_dir = "S$fix"."_exp$hr->{experiment}";
        if ($fwd_file && $rev_file) {
            my $crispresso_cmd;
            if ($single) {
                my $file;
                if ($hr->{strand} eq '-') {
                    $file = $rev_file;
                }
                else {
                    $file = $fwd_file;
                }
                $crispresso_cmd = "$crispresso -w 50 --hide_mutations_outside_window_NHEJ --save_also_png "
                        ." -o $output_dir"
                        ." -r1 $file"
                        ." -a $hr->{amplicon}"
                        ." -g $crispr_site";
            } else {
                $crispresso_cmd = "$crispresso -w 50 --hide_mutations_outside_window_NHEJ --save_also_png "
                        ." -o $output_dir"
                        ." -r1 $fwd_file"
                        ." -r2 $rev_file"
                        ." -a $hr->{amplicon}"
                        ." -g $crispr_site";
            }

            if ($hr->{hdr}) {
                $crispresso_cmd = $crispresso_cmd . " -e $hr->{hdr}";
            }

            my $bsub_cmd = 'bsub -n1 -q normal -G team87-grp -M2000 -R"select[mem>2000] rusage[mem=2000] span[hosts=1]"'
                    ." -eo $output_dir/job.err"
                    ." -oo $output_dir/job.out"
                    ." $crispresso_cmd";

            say "Running command $bsub_cmd";
            my $cmd_output = `$bsub_cmd`;
            say "Command output: $cmd_output";
        } else {
            say "No file found for barcode $barcode";
        }
    }
}

