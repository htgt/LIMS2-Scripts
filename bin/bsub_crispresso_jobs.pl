#! /usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;
use Path::Class;
use Getopt::Long;

# Inputs:
#   samples: csv file containing 5 columns: Experiment ID, crispr seq, crispr strand, amplicon seq, barcode range
#   dir    : directory containing input fastq files (1 per barcode with names containing _S<barcode number>_ and R1 or R2)

# NB: Experiment ID is just a unique identifier for the combination of amplicon and crispr used.
# It is not a LIMS2 experiment ID at this stage in development.

my $crispresso = $ENV{CRISPRESSO_CMD} || "/nfs/team87/farm3_lims2_vms/software/python_local/bin/CRISPResso";

GetOptions(
    'samples=s' => \my $samples_file_name,
    'dir=s'      => \my $dir_name,
);

# Ouputs:
# 1 directory per barcode/experiment combination
# directory names in format S<barcode number> _expt<experiment ID>

# bsub -q normal -G team87-grp -eo S1_expA_farm/job.err -oo S1_expA_farm/job.out -M500 -R"select[mem>500] rusage[mem=500]" ./bin/CRISPResso -w 30 --hide_mutations_outside_window_NHEJ --ignore_substitutions -o S1_expA -r1 real_data/orig_fastq/Homo-sapiens_S1_L001_R1_001.fastq -a GAAAGTCCGGAAAGACAAGGAAGgaacacctccacttacaaaagaagataagacagttgtcagacaaagccctcgaaggattaagccagttaggattattccttcttcaaaaaggacagatgcaaccattgctaagcaactcttacagag -g CTCGAAGGATTAAGCCAGTT
my $dir = dir($dir_name);
my @files = $dir->children;
my @samples_info = file($samples_file_name)->slurp(chomp => 1);
foreach my $line (@samples_info){
	my ($exp_id, $gene, $crispr_seq, $crispr_strand, $amplicon, $barcode_range) = split /\s*,\s*/, $line;
    my ($start,$end) = split /\s*-\s*/, $barcode_range;

    # remove PAM, revcom if crispr site is on negative strand
    my $crispr_site;
    #if($crispr_strand eq "+"){
    $crispr_site = substr($crispr_seq,0,20);
    #}
    #elsif($crispr_strand eq "-"){
    #	$crispr_seq = reverse scalar $crispr_seq;
    #	$crispr_seq =~ tr/ATCG/TAGC/;
    #	$crispr_site = substr($crispr_seq,0,20);
    #}

    my $amplicon_start = substr($amplicon,0,150);

    my @barcodes = $start..$end;
    foreach my $barcode (@barcodes){
    	my $bc_in_name = "_S".$barcode."_";
        my $file;
        if($crispr_strand eq "-"){
            ($file) = grep { $_=~ /$bc_in_name.*R2/ } @files;
        } else {
            ($file) = grep { $_=~ /$bc_in_name.*R1/ } @files;
        }
        if($file){
        	my $output_dir = "S$barcode"."_exp$exp_id";
            my $crispresso_cmd = "$crispresso -w 30 --hide_mutations_outside_window_NHEJ --ignore_substitutions"
                                 ." -o $output_dir"
                                 ." -r1 $file"
                                 ." -a $amplicon_start"
                                 ." -g $crispr_site";

            my $bsub_cmd = 'bsub -n1 -q normal -G team87-grp -M1000 -R"select[mem>1000] rusage[mem=1000] span[hosts=1]"'
                           ." -eo $output_dir/job.err"
                           ." -oo $output_dir/job.out"
                           ." $crispresso_cmd";

            say "Running command $bsub_cmd";
            my $cmd_output = `$bsub_cmd`;
            say "Command output: $cmd_output";
        }
        else{
        	say "No file found for barcode $barcode";
        }
    }
}

