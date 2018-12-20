#!/usr/local/bin/bash
PYTHONPATH=/nfs/team87/farm3_lims2_vms/software/python_local/lib/python/
PATH=/nfs/team87/farm3_lims2_vms/software/crispresso_dependencies/bin:$PATH
crispresso=${CRISPRESSO_CMD:-/nfs/team87/farm3_lims2_vms/software/python_local/bin/CRISPResso}
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
fwd=`find . -maxdepth 1 -name "*_S${LSB_JOBINDEX}_*R1*"`
rev=`find . -maxdepth 1 -name "*_S${LSB_JOBINDEX}_*R2*"`
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

let fix=$LSB_JOBINDEX-$offset
JOB="$crispresso -w 50 --hide_mutations_outside_window_NHEJ --save_also_png -o S${fix}_exp$name -a $amplicon -g $crispr" 
if [ $single -ne 0 ]; then
    if [ $reverse -ne 0 ]; then
        JOB="$JOB -r1 $rev"
    else
        JOB="$JOB -r1 $fwd"
    fi
else
    JOB="$JOB -r1 $fwd -r2 $rev"
fi
if [ $hdr ]; then
    JOB="$JOB -e $hdr"
fi
echo $LSB_JOBINDEX
echo $JOB
$JOB
