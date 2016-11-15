#!/usr/bin/env perl

use strict;
use warnings;

use LIMS2::Model;
use Log::Log4perl qw( :easy );
use Data::Dumper;
use feature qw(say);
use Try::Tiny;

my $output_file = $ARGV[0] or die "Please provide the output file name as an argument to the script";
open (my $fh, ">", $output_file) or die "Cannot open file $output_file for writing";

my $model = LIMS2::Model->new( user => 'lims2' );

my $species = 'Human';

my @plates = $model->schema->resultset('Plate')->search({
	type_id => 'EP_PICK',
	species_id    => $species,
	})->all;

my @headings = qw(
    plate_id
    plate_name
    well_id
    well_name
    qc_well_id
    qc_run_id
    qc_run_date
    qc_well_accepted
    crispr_damage_type
    variant_size
);

print $fh join "\t", @headings;
print $fh "\n";

foreach my $plate (@plates){
	say "Getting QC for plate $plate";
	foreach my $well ($plate->wells){
		my $has_qc = 0;
		foreach my $qc_well ($well->crispr_es_qc_wells){
			$has_qc = 1;
            my @data = (
                $plate->id,
                $plate->name,
                $well->id,
                $well->name,
                $qc_well->id,
                $qc_well->crispr_es_qc_run->id,
                $qc_well->crispr_es_qc_run->created_at->iso8601,
                $qc_well->accepted,
                ( $qc_well->crispr_damage_type_id // '' ),
                ( $qc_well->variant_size // '' ),
            );
            print $fh join "\t", @data;
            print $fh "\n";
		}

		unless($has_qc){
			my @data = (
                $plate->id,
                $plate->name,
                $well->id,
                $well->name,
                "no qc",
			);
			print $fh join "\t", @data;
			print $fh "\n";
		}
	}
}
close $fh;
