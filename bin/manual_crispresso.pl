#! /usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;
use Path::Class;
use Getopt::Long;

# Inputs:
#   samples: csv file containing 6 columns: Experiment ID, crispr seq, crispr strand, amplicon seq, min index, max index
#   dir    : directory containing input fastq files (1 per barcode with names containing _S<barcode number>_ and R1 or R2)

# NB: Experiment ID is just a unique identifier for the combination of amplicon and crispr used.
# It is not a LIMS2 experiment ID at this stage in development.

my $crispresso = $ENV{CRISPRESSO_CMD} || "/nfs/team87/farm5/software/python_local/bin/CRISPResso";

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

my $dir = dir($dir_name);
my @files = $dir->children;
my @samples_info = file($samples_file_name)->slurp(chomp => 1);
OUTER: foreach my $line (@samples_info) {
	my ($exp_id, $gene, $crispr_seq, $crispr_strand, $amplicon, $min_index, $max_index, $hdr) = split /\s*,\s*/, $line;

    if ($min_index eq 'min_index') {next OUTER};

    my @barcodes = $min_index .. $max_index;
    my $crispr_site = substr($crispr_seq,0,20);

    foreach my $barcode (@barcodes) {
    	my $bc_in_name = "_S".$barcode."_";
        my ($fwd_file) = grep { $_=~ /$bc_in_name.*R1/ } @files;
        my ($rev_file) = grep { $_=~ /$bc_in_name.*R2/ } @files;

        my $fix = $barcode - $offset;
        my $output_dir = "S$fix"."_exp$exp_id";

        if ($fwd_file && $rev_file) {
            my $crispresso_cmd;
            if ($single) {
                my $file;
                if ($crispr_strand eq '-') {
                    $file = $rev_file;
                } else {
                    $file = $fwd_file;
                }

                $crispresso_cmd = "$crispresso -w 50 --hide_mutations_outside_window_NHEJ --save_also_png "
                     ." -o $output_dir"
                     ." -r1 $file"
                     ." -a $amplicon"
                     ." -g $crispr_site";
            } else {
                $crispresso_cmd = "$crispresso -w 50 --quality --hide_mutations_outside_window_NHEJ --save_also_png "
                     ." -o $output_dir"
                     ." -r1 $fwd_file"
                     ." -r2 $rev_file"
                     ." -a $amplicon"
                     ." -g $crispr_site";
            }

            if ($hdr) {
                $crispresso_cmd = $crispresso_cmd . " -e $hdr";
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

