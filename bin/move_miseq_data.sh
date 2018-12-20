#!/usr/local/bin/bash
while getopts "p:r:" opt; do
    case "$opt" in
        p) processed=$OPTARG ;;
        r) raw=$OPTARG ;;
    esac
done
die() { echo "$*" 1>&2; exit 1; }
if [[ ! $processed ]] || [[ ! $raw ]]; then
    die "Usage: $0 -p <processed destination> -r <raws destination>"
fi
if [[ ! -d $processed ]]; then
    die "Processed destination $processed doesn't exist"
fi
if [[ ! -d $raw ]]; then
    die "Raw destination $raw doesn't exist"
fi

echo "Storing raw files in $raw"
find . -maxdepth 1 -name "*.fastq.gz" -exec mv -t "$raw" {} +
echo "Deleting compressed files"
find ./S*_exp*/C*/ -name \*.gz -delete
echo "Storing processed files in $processed"
find . -maxdepth 1 -type d -name "S*" -exec mv -t "$processed" {} +
cp summary.csv "$processed/summary.csv" 
echo "Setting write permissions"
find "$processed" -type d -exec chmod 775 {} \;
chmod 664 "$processed/summary.csv"
