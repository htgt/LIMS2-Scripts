#!/usr/bin/env perl

use strict;
use warnings;

use Log::Log4perl qw( :easy );
use Data::Dumper;
use Data::Compare;
use feature qw(say);
use LIMS2::Model;
use Try::Tiny;
use Text::CSV;

my $model = LIMS2::Model->new( user => 'lims2' );

#Tried with array. Bad idea. Hogs memory and processing blocks
#Seperated into two queries to lower the processing load
#Distinct will ignore repeating data
my $plate_rs_ref = $model->schema->resultset('Summary')->search(
    {
        design_species_id => 'Mouse',
        ep_pick_plate_name => {'!=',undef}
    },
    {
        columns => [qw/
            design_id
        /],
        distinct => 1,
    } );

my $design_info;
my @designed_plates;
my $summary;
my $count = $plate_rs_ref->count;

#CSV printing
my $csv = Text::CSV->new ( { binary => 1 , auto_diag => 1, eol => "\n"} )
                 or die "Cannot use CSV: ".Text::CSV->error_diag ();
#Insert column names to top of csv
$csv->column_names('EPD_Name','EPD Parent_Cell_Line', 'EPD_Backbone', 'EPD_Cassette_Name', 'EPD_Cassette_Promoter','Vector_Name', 'Vector_Backbone',
        'Gene_Marker_Symbol', 'Design_ID', 'Design_Type','Loxp_Start','Loxp_End','Cassette_Start', 'Cassette_End', 'Homology_Start', 'Homology_End');
$csv->bind_columns();

#Opens the file and creates a file handle
my $fh;
open $fh, ">", "results.csv" or die "Failed: $!";
$csv->print ($fh, [$csv->column_names]);
close $fh;

#Iterate through matches
while ( my $plate_rs = $plate_rs_ref->next ) {
    #Finds the corresponding summary to retrieve all necessary data
    #Column specification reduces the number of results with distinct 
    $summary = $model->schema->resultset('Summary')->search(
        {
            design_id => $plate_rs->{_column_data}->{design_id},
            design_species_id => 'Mouse',
            ep_pick_plate_name => {'!=',undef}
        },
        {
            columns => [qw/
                ep_pick_plate_name
                ep_pick_well_name
                ep_first_cell_line_name
                int_backbone_name
                final_cassette_name
                final_cassette_promoter
                final_plate_name
                final_backbone_name
                design_gene_symbol
                design_id
                design_type
            /],
            distinct => 1,
        }
    );
    #For each summary relating to the design id, find the design info and bring to form the csv 
    while ( my $summary_result = $summary->next ) {
        if (exists $summary_result->{_column_data}){
            $design_info = retrieve_design_info($model, $plate_rs);
            @designed_plates = build_plate_array($summary_result, $design_info, @designed_plates);
            print_to_csv($csv, @designed_plates);
        }
    }
    print $count, "\n";
    $count = $count - 1;

}


#Iterates through printing each line to the csv
sub print_to_csv {
    my ($csv_f, @plates) = @_;
    pop(@plates);
    # >> signals append to file
    open my $file_handle, ">>", "results.csv" or die "Failed: $!";
    $csv_f->print($file_handle, \@plates);
    close $file_handle;

    return;
}

#Creates a design templete. Using that template query for relating design info
sub retrieve_design_info{
    my ($model_db, $plate) = @_;
    my $design_rs = $model_db->schema->resultset('Design')->find($plate->{_column_data}->{design_id});
    my $design_info_rs = LIMS2::Model::Util::DesignInfo->new( design => $design_rs );
    return $design_info_rs;
}

#Places all the data into the array for ease of printing
sub build_plate_array {
    my ($plate, $design_info_ref, @plates_info) = @_;
    my @plate_array;
    my $promoter_state;

    my $epd_name = join('_', $plate->{_column_data}->{ep_pick_plate_name}, $plate->{_column_data}->{ep_pick_well_name});

#If undefined in the hash ref, set state to not defined else decide promoter state
#Exists will count an undef, defined won't
    if (defined $plate->{_column_data}->{final_cassette_promoter}){
        if ( $plate->{_column_data}->{final_cassette_promoter} == 1 ){
            $promoter_state = "Promoter driven";
        }
        else{
            $promoter_state = "Promoterless";
        }
    }
    else{
        $promoter_state = "Not defined";
    }

    my $loxp_start;
    my $loxp_end;

    if ($design_info_ref->_build_loxp_start && $design_info_ref->_build_loxp_end){
        $loxp_start = $design_info_ref->_build_loxp_start;
        $loxp_end = $design_info_ref->_build_loxp_end;
    }
    else{
        $loxp_start = $plate->{_column_data}->{design_type};
        $loxp_end = $plate->{_column_data}->{design_type};
    }

    @plate_array = (
        $epd_name,
        $plate->{_column_data}->{ep_first_cell_line_name},
        $plate->{_column_data}->{int_backbone_name},
        $plate->{_column_data}->{final_cassette_name},
        $promoter_state,
        $plate->{_column_data}->{final_plate_name},
        $plate->{_column_data}->{final_backbone_name},
        $plate->{_column_data}->{design_gene_symbol},
        $plate->{_column_data}->{design_id},
        $plate->{_column_data}->{design_type},
        $loxp_start,
        $loxp_end,
        $design_info_ref->_build_cassette_start,
        $design_info_ref->_build_cassette_end,
        $design_info_ref->_build_homology_arm_start,
        $design_info_ref->_build_homology_arm_end,
    );
    return @plate_array;
}
