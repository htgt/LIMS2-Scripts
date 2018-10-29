#!/usr/local/bin/bash
DEST=$1
echo "Destination: '$DEST'"
mv S* $DEST 
cp summary.csv $DEST/summary.csv 
