#!/usr/local/bin/bash
PERL5LIB=/nfs/team87/farm3_lims2_vms/software/perl/lib/perl5
export BASESPACE_TOKEN=$1
echo "Job: $LSB_JOBINDEX"
SAMPLE="$(head -$LSB_JOBINDEX samples.txt | tail -1)"
echo "Sample: $SAMPLE"
BIN="$(readlink -f "$(dirname "$0")")"
$BIN/download_basespace.pl --sample=$SAMPLE --path=.
