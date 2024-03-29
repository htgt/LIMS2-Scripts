#!/usr/local/bin/bash
export PYTHONPATH=/nfs/team87/farm5/software/python_local/lib/python/
export PATH=/nfs/team87/farm5/software/crispresso_dependencies/bin:$PATH
export crispresso=${CRISPRESSO_CMD:-/nfs/team87/farm5/software/python_local/bin/CRISPResso}
single=0
reverse=0
offset=0
end() { echo "$*" 1>&2 ; exit 0; }
die()  { echo "$*" 1>&2 ; exit 1; }
while getopts "rso:a:g:n:e:" opt; do
    case "$opt" in
    s) single=1 ;;
    r) reverse=1 ;;
    o) offset=$OPTARG ;;
    a) amplicon=$OPTARG ;;
    g) crispr=$OPTARG ;;
    n) name=$OPTARG ;;
    e) hdr=$OPTARG ;;
    esac
done
let illumina_index=$LSB_JOBINDEX+$offset
fwd=`find . -maxdepth 1 -name "*_S${illumina_index}_*R1*"`
rev=`find . -maxdepth 1 -name "*_S${illumina_index}_*R2*"`
if [ ! $amplicon ]; then
    die "No amplicon specified"
fi
if [ ! $crispr ]; then
    die "No CRISPR guide specified"
fi
if [ ! $fwd ]; then
    end "Could not find forward file"
fi
if [ ! $rev ]; then
    end "Could not find reverse file"
fi

JOB="$crispresso -w 50 --quality --hide_mutations_outside_window_NHEJ --save_also_png -o S${LSB_JOBINDEX}_exp$name -a $amplicon -g $crispr"
if [ $single -ne 0 ]; then
    if [ $reverse -ne 0 ]; then
        JOB="$JOB -r1 $rev"
    else
        JOB="$JOB -r1 $fwd"
    fi
else
    JOB="$JOB -r1 $fwd -r2 $rev"
fi
if [ ! -z $hdr ]; then
    JOB="$JOB -e $hdr"
fi
echo $LSB_JOBINDEX
echo $JOB
$JOB
