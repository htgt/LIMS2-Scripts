This folder contains a work in progress for refactoring of scripts for 
the crispr conditional strategy: LIMS2-Scripts/crispr_conditional/
step1_target_finder.pl
step2_target_primers.pl
step3_crispr_selector.pl

step1_target_finder.pl locates the insertion sites for a list of genes. 
    INPUT: gene_list_file
    OUTPUT: targets_file
step2_target_primers.pl finds primers for each insertion site.
    INPUT: targets_file
    OUTPUT: targets_primers_file
step3_crispr_selector.pl finds and ranks crisprs for each target.
    INPUT: targets_primers_file
    OUTPUT: targets_primers_crisprs_file

The refactoring still needs work to be done, but is functional and it is recommended to use it in new computations.
For recomputations I recommend using the original scrips to make sure nothing has changed. 

Original scripts are in: LIMS2-Scripts/bin/
bk_target_finder.pl
bk_target_crispr_finder.pl
crispr_conditional_primers.pl
crispr_selector.pl

bk_target_finder.pl locates the insertion sites for a list of genes.
    INPUT: gene_list_file
    OUTPUT: targets_file
bk_target_crispr_finder.pl finds crisprs for each insertion site.
    INPUT: targets_file
    OUTPUT: targets_crisprs_file
crispr_conditional_primers.pl finds primers for each insertion site.
    INPUT: targets_file
    OUTPUT: targets_primers_file
crispr_selector.pl finds and ranks crisprs for each target.
    INPUT: targets_primers_file
    OUTPUT: targets_primers_crisprs_file 

A list of genes is required (all_human_genes.csv and all_mouse_genes.csv)
together with primer3 config (primer3_crispr_conditional_la_config.yaml and primer3_crispr_conditional_ra_config.yaml)
These can be modified to meet primer requirements.

Data files of runs have been archived and can be found at
/warehouse/team229_wh01/tg6/crispr_conditional/Backup

