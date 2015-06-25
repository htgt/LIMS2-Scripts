#!/usr/bin/env perl

use strict;
use warnings;

use LIMS2::Model;
use Log::Log4perl qw( :easy );
use Data::Dumper;
use feature qw(say);
use Try::Tiny;

=head

Quick n dirty script to compare the crispr damage calls we have made in LIMS2 to those
made by VEP as our calls did not take into account intronic vs exonic damage

Output is tab delimited text file giving QC well info, our call and VEP call for the
gene of interest

Human only

=cut

my $model = LIMS2::Model->new( user => 'lims2' );

#my @qc_wells = $model->schema->resultset('CrisprEsQcWell')->search({ id => 9465 })->all;
my @qc_wells = $model->schema->resultset('CrisprEsQcWell')->all;

my $gene_finder = sub{ $model->find_genes( @_ ) };

my @display_items = qw(es_qc_well_id gene crispr_id is_crispr_pair is_crispr_group damage_type);
my @headings = qw(well_name well_accepted es_qc_well_accepted);
push @headings, @display_items, "vep_call", "qc_run_identifier";
my @lines;

say join "\t", @headings;

foreach my $qc_well (@qc_wells){
	next unless $qc_well->crispr_es_qc_run->species_id eq 'Human';

	my $data;
	try { $data = $qc_well->format_well_data($gene_finder) };

	unless ($data){
		print STDERR "Could not get well data for ".$qc_well->well->as_string."\n";
		say $qc_well->well->as_string;
		next;
	}

    my @values = (
        $qc_well->well->as_string,
        $qc_well->well->is_accepted,
        $qc_well->accepted,
    );

    push @values, @{ $data }{@display_items};

    my $vep_call = get_vep_call($data->{vep_output}, $data->{gene_ids}, $data->{gene});
    push @values, $vep_call, $qc_well->crispr_es_qc_run_id;

    say join "\t", @values;
}

say 'all crispr QC wells examined';

sub get_vep_call{
	my ($vep_output, $gene_ids, $gene_symbol) = @_;

	return 'no vep output' unless $vep_output;

	my ($gene_id) = @$gene_ids;

	my $gene_variant;
    if($gene_id){
    	($gene_variant) = grep { $_ =~ /HGNC_ID=$gene_id/ } split("\n", $vep_output);
    }

    # Older VEP output format (seen in v75)
    unless($gene_variant){
    	if($gene_symbol){
    	    ($gene_variant) = grep { $_ =~ /SYMBOL=$gene_symbol/ } split("\n", $vep_output);
        }
    }

	return "no vep output for $gene_id/$gene_symbol" unless $gene_variant;

    # Consequence is at index 6 (FIXME: do not hard code this in production version)
	my @items = split ("\t",$gene_variant);
	return $items[6];
}