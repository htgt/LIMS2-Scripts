#!/usr/bin/env perl

use strict;
use warnings;

use LIMS2::Model;
use Log::Log4perl qw( :easy );
use Data::Dumper;
use feature qw(say);

open (my $output, ">", "intronic_crisprs.csv");
print $output join ",",qw(crispr_id crispr_chr crispr_start crispr_end crispr_wells assembly_wells ep_wells gene_name);
print $output "\n";
my $model = LIMS2::Model->new( user => 'lims2' );

# I think we only care about crisprs we have started work on, e.g. plated crisprs
my $crispr_rs = $model->schema->resultset('ProcessCrispr')->search_related('crispr', { crispr_loci_type_id => 'Intronic', species_id => 'Human' });
my @crisprs = $crispr_rs->search({},{ distinct => 1, order_by => { '-asc' => 'id' } });

say "We have ".scalar(@crisprs)." crisprs to check";
my $gene_finder = sub { $model->find_genes( @_ ); };
my $counter=0;
foreach my $crispr (@crisprs){
    $counter++;

    say $counter;

    my @crispr_wells;
    foreach my $process_crispr ( $crispr->process_crisprs ){
    	my $process = $process_crispr->process;
    	push @crispr_wells, $process->output_wells;
    }
    my $well_names = join "/",( map { $_->as_string } @crispr_wells );

    my @assembly_wells = map { $_->descendants_of_type('ASSEMBLY') } @crispr_wells;
    my $assembly_well_types = {};
    my $gene_name;
    foreach my $well (@assembly_wells){
        my ($assembly_process) = $well->parent_processes;
        $assembly_well_types->{ $well->as_string } = $assembly_process->type_id;
        unless($gene_name){
            my $design = $well->design;
            ($gene_name) = $design->gene_symbols($gene_finder);
        }
    }
    my $assembly_names = join "/", (map { $_->as_string."(".$assembly_well_types->{$_->as_string}.")" } @assembly_wells );
    my @ep_wells = map { $_->descendants_of_type('CRISPR_EP') } @crispr_wells;
    my $ep_names = join "/", (map { $_->as_string } @ep_wells );

    print $output join ",",
    	$crispr->id,
        $crispr->chr_name,
    	$crispr->start,
    	$crispr->end,
    	$well_names,
        $assembly_names,
        $ep_names,
        $gene_name;

    print $output "\n";
}

