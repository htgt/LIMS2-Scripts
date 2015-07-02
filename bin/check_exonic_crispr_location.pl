#!/usr/bin/env perl

use strict;
use warnings;

use LIMS2::Model;
use Log::Log4perl qw( :easy );
use Data::Dumper;
use feature qw(say);

open (my $no_exon_list, ">", "no_exon_found.csv");
open (my $output, ">", "boundary_crisprs.csv");
print $output join ",",qw(crispr_id crispr_start crispr_end crispr_wells assembly_wells ep_wells gene_name exon_id exon_start exon_end);
print $output "\n";
my $model = LIMS2::Model->new( user => 'lims2' );

# I think we only care about crisprs we have started work on, e.g. plated crisprs
my @crisprs = $model->schema->resultset('ProcessCrispr')->search_related('crispr', { crispr_loci_type_id => 'Exonic', species_id => 'Human' })->all;
#my @crisprs = $model->schema->resultset('Crispr')->search({ id => { -in => [191466,192182] } })->all;

say "We have ".scalar(@crisprs)." crisprs to check";
my $counter=0;
foreach my $crispr (@crisprs){
    $counter++;
 	my $slice = $model->ensembl_slice_adaptor('human')->fetch_by_region(
		'chromosome',
		$crispr->chr_name,
		$crispr->start,
		$crispr->end,
	);

    my ($exon) = @{ $slice->get_all_Exons };
    
    # Check we found the exon
    unless($exon){
    	warn "$counter: No exon found for crispr ".$crispr->id;
    	print $no_exon_list $crispr->id;
    	print $no_exon_list "\n";
    	next;
    }

    # Ensure exon start is always before end (not sure how ensembl coords work)
    my ($exon_start,$exon_end) = sort ($exon->seq_region_start, $exon->seq_region_end);
    if($crispr->start < $exon_start or $crispr->end > $exon_end){
    	say "$counter: Crispr ".$crispr->id." overlaps intron exon boundary";
    	my @crispr_wells;
    	foreach my $process_crispr ( $crispr->process_crisprs ){
    		my $process = $process_crispr->process;
    		push @crispr_wells, $process->output_wells; 
    	}
    	my $well_names = join "/",( map { $_->as_string } @crispr_wells );

    	my @assembly_wells = map { $_->descendants_of_type('ASSEMBLY') } @crispr_wells;
        my $assembly_well_types = {};
        foreach my $well (@assembly_wells){
        	my ($assembly_process) = $well->parent_processes;
        	$assembly_well_types->{ $well->as_string } = $assembly_process->type_id;
        }
    	my $assembly_names = join "/", (map { $_->as_string."(".$assembly_well_types->{$_->as_string}.")" } @assembly_wells );
    	my @ep_wells = map { $_->descendants_of_type('CRISPR_EP') } @crispr_wells;
    	my $ep_names = join "/", (map { $_->as_string } @ep_wells );

    	print $output join ",", 
    	$crispr->id, 
    	$crispr->start, 
    	$crispr->end, 
    	$well_names,
        $assembly_names,
        $ep_names,
        $exon->get_nearest_Gene->external_name,
    	$exon->stable_id, 
    	$exon_start, 
    	$exon_end;

        print $output "\n";
    }
    else{
    	say "$counter: Crispr ".$crispr->id." is fully exonic";
    }
}
