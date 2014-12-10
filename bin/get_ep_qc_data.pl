#!/usr/bin/env perl

use strict;
use warnings;

use feature qw( say );
use Data::Dumper;

use LIMS2::Model;
use Text::CSV;

my $model = LIMS2::Model->new( user => 'lims2' );

open my $fh, "<", $ARGV[0] or die "$!";
open my $out_fh, ">", "crispr_ep_data.csv" or die "$!";

my $csv = Text::CSV->new( { eol => "\n" });
$csv->column_names( $csv->getline( $fh ) );

my @columns = $csv->column_names;
push @columns, (
    "Num GR bands",
    "Num GF bands",
    "Num TR bands",
    "Num valid LR reads",
    "Num valid R1R/R2R reads"
);

$csv->print( $out_fh, \@columns );

while ( my $data = $csv->getline_hr( $fh )) {
    my $plate_name = $data->{"Crispr EP Plate Name"};
    my $well_name  = $data->{"Crispr EP Well Name"};

    #my $well = $model->retrieve_well( { plate_name => "HEP00994", well_name => "A02" } );
    my $well = $model->retrieve_well( { plate_name => $plate_name, well_name => $well_name } );

    my @epds = get_epd_wells( $well );

    #say $_->plate->name for @epds;

    my %summary_data;
    for my $w ( @epds ) {
        my $qc_result = $w->well_qc_sequencing_result;

        if ( $qc_result ) {
            #say Dumper( $qc_result->as_hash );
            my $qc_data = $qc_result->as_hash;
            my @valid_primers = split ",", $qc_data->{valid_primers};
            for my $p ( @valid_primers ) {
                $summary_data{num_gr}++ if $p =~ /^GR/;
                $summary_data{num_gf}++ if $p =~ /^GF/;
                $summary_data{num_tr}++ if $p =~ /^TR/;
                $summary_data{num_lr}++ if $p =~/^LR$/;
                $summary_data{num_r}++ if $p =~ /R1R|R2R/;
            }
        }
        else {
            $summary_data{num_no_seq}++;
            #say $w->id . " " . $w->name . " " . $w->plate->name;
            for my $band ( $w->well_primer_bands ) {
                next unless $band->pass;
                my $id = $band->primer_band_type->id;
                $summary_data{num_gr}++ if $id =~ /^gr/;
                $summary_data{num_gf}++ if $id =~ /^gf/;
                $summary_data{num_tr}++ if $id =~ /tr_pcr/;
            }

            #die "Well has no QC data!";
        }
    }

    my @output_data = @{ $data }{ $csv->column_names };
    push @output_data, (
        $summary_data{num_gr} || 0,
        $summary_data{num_gf} || 0,
        $summary_data{num_tr} || 0,
        $summary_data{num_lr} || 0,
        $summary_data{num_r} || 0,
    );

    $csv->print( $out_fh, \@output_data );

    say Dumper( \%summary_data );
}

sub get_epd_wells {
    my ( $well ) = @_;

    my @pick_wells;
    my $descendants = $well->descendants->depth_first_traversal( $well, 'out' );
    while( my $d = $descendants->next ) {
        if ( $d->plate->type_id eq 'EP_PICK' ) {
            push @pick_wells, $d;
        }
    }

    say "Found " . scalar( @pick_wells ) . " EPD wells";

    return @pick_wells;
}