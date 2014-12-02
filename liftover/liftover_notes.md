# Steps
- Run following script to download current coordinates into bed files:
    - specify species and current assembly
    - `fetch_coords.pl --species Human --assembly GRCh37`
- Multiple bed files will be created with coordinates from:
    - `design_oligo_loci`
    - `genotyping_primer_loci`
    - `crispr_loci`
    - `crispr_primer_loci`
- Run liftover on the 4 bed files ( See section below ).
- Run the `add_new_loci.pl` script for each of the 4 tables
    - must specify species, new assembly, data-file and table resultset name.
    - run initially without commit flag, if everything goes smoothly run script with commit flag.
    - `add_new_loci.pl --species Human --assembly GRCh38 --resultset DesignOligoLocus --data-file DesignOligoLocus.new.bed --commit`

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
