#!/usr/bin/env perl

use strict;
use warnings;

use LIMS2::Model;
use Log::Log4perl qw( :easy );
use Data::Dumper;
use feature qw(say);

my $species = 'Mouse';
my $external_db = "MGI";

open (my $no_exon_list, ">", $species."_no_exon_found.csv");
open (my $output, ">", $species."_boundary_crisprs.csv");
print $output join ",",qw(crispr_id lims2_crispr_type chromosome crispr_start crispr_end crispr_wells assembly_wells ep_wells gene_name exon_id exon_start exon_end sponsors);
print $output "\n";
my $model = LIMS2::Model->new( user => 'lims2' );
my $gene_finder = sub { $model->find_genes( @_ ); };

# I think we only care about crisprs we have started work on, e.g. plated crisprs
my @crisprs = $model->schema->resultset('ProcessCrispr')->search_related('crispr', { species_id => $species })->all;
#my @crisprs = $model->schema->resultset('Crispr')->search({ id => { -in => [191466,192182] } })->all;

say "We have ".scalar(@crisprs)." crisprs to check";
my $counter=0;
foreach my $crispr (@crisprs){
    $counter++;
 	my $slice = $model->ensembl_slice_adaptor(lc($species))->fetch_by_region(
		'chromosome',
		$crispr->chr_name,
		$crispr->start,
		$crispr->end,
	);

    my ($gene) = @{ $slice->get_all_Genes };

    unless($gene){
        say "$counter: no gene found for crispr ".$crispr->id;
        my ($well_names, $assembly_names, $ep_names) = get_well_names($crispr);

        print $output join ",",
        $crispr->id,
        $crispr->crispr_loci_type_id,
        $crispr->chr_name,
        $crispr->start,
        $crispr->end,
        $well_names,
        $assembly_names,
        $ep_names,
        'no gene',
        'no gene',
        '',
        '',
        '';

        print $output "\n";
        next;
    }

    my %canonical_exons = map { $_->stable_id => 1 } @{ $gene->canonical_transcript->get_all_Exons };

    my ($exon) = grep { exists $canonical_exons{ $_->stable_id } } @{ $slice->get_all_Exons };

    # Check we found the exon
    if($exon){
        # Ensure exon start is always before end (not sure how ensembl coords work)
        my ($exon_start,$exon_end) = sort ($exon->seq_region_start, $exon->seq_region_end);
            my ($well_names, $assembly_names, $ep_names) = get_well_names($crispr);
            my $sponsors  = get_sponsors($gene);

            print $output join ",",
            $crispr->id,
            $crispr->crispr_loci_type_id,
            $crispr->chr_name,
            $crispr->start,
            $crispr->end,
            $well_names,
            $assembly_names,
            $ep_names,
            $gene->external_name,
            $exon->stable_id,
            $exon_start,
            $exon_end,
            $sponsors;

        if($crispr->start < $exon_start or $crispr->end > $exon_end){
        	say "$counter: Crispr ".$crispr->id." overlaps intron exon boundary";

            print $output ",boundary overlap\n";
        }
        else{
        	say "$counter: Crispr ".$crispr->id." is fully exonic (".$exon->stable_id.")";

            print $output ",fully exonic\n";
        }
    }
    else{
        my ($well_names, $assembly_names, $ep_names) = get_well_names($crispr);
        my $sponsors = get_sponsors($gene);
        # No exon found. this is not good!
        say "$counter: Crispr ".$crispr->id." does not overlap a coding exon in canonical transcript";

        print $output join ",",
        $crispr->id,
        $crispr->crispr_loci_type_id,
        $crispr->chr_name,
        $crispr->start,
        $crispr->end,
        $well_names,
        $assembly_names,
        $ep_names,
        $gene->external_name,
        'no exon found',
        '',
        '',
        $sponsors;

        print $output "\n";
    }
}

sub get_well_names {
    my $crispr = shift;

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

    return ($well_names, $assembly_names, $ep_names);
}

sub get_sponsors{
    my $gene = shift;
    my $gene_info = $model->find_gene( { search_term => $gene->external_name, species => $species } );
    if($gene_info){
        my $gene_id = $gene_info->{gene_id};
        my $projects_rs = $model->schema->resultset('Project')->search(
            {
                gene_id => $gene_id,
            }
            ,
            {
                order_by => 'id',
                join => 'project_sponsors',
                distinct => 1,
            }
        );
        my %sponsors = map { $_ => 1 } map { $_->sponsor_ids } $projects_rs->all;
        my $sponsor_string = join "/", sort keys %sponsors;
        return $sponsor_string;
    }
    return "";
}