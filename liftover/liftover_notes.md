# UCSC Liftover Tool
Work out new assembly coordinates given the old assembly coordinates.
Download from: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/

Grab chain file for our conversion from: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/
 hg19ToHg38.over.chain.gz

## Liftover Tool Steps
- Grab current coordinates for species we are lifting over
- Create a bed file of this data, the id of each oligo should be something we can use to insert data back into LIMS2
    - e.g. design oligo id
- Grab correct chain file from UCSC, need to be for the right species and the right transition
    - e.g. mm9ToMm10.over.chain
- run liftover tool on bed file with the chain file
    - `liftOver input.bed chain_file output.bed unmatched.bed`
- check the unmatched.bed file for any liftover failures ( fix if you can )
- Insert the new oligo coordinates from the output.bed file into LIMS2
