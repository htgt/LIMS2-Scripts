#!/usr/bin/env perl

use strict;
use warnings;

use File::Slurp;
use Text::CSV;
use LIMS2::Model;

my $model = LIMS2::Model->new ( user => 'lims2' );

my @output_data; 
my @headers = ('Barcode', 'CRISPR seq', 'PAM right', 'Chromosome');
push @output_data, \@headers;
my $file = $ARGV[0];
my @input_data = read_file($file, {chomp => 1});
foreach my $barcode (@input_data) {
    # remove any whitespace characters
    ($barcode) = $barcode =~ m/([A-Z0-9]+)/;
    my $well = $model->schema->resultset('Well')->find( { barcode => $barcode } );
    if (! $well) {
        print "Could not find well with barcode: $barcode!\n";
        my @output_row = ($barcode, '', '');
        push @output_data, \@output_row;
        next;
    }
    foreach my $crispr_well (map {$_->crispr->as_hash} $well->parent_crispr_wells) {
        my @output_row = ($barcode, $crispr_well->{fwd_seq}, $crispr_well->{pam_right}, "$crispr_well->{locus}->{chr_name}:$crispr_well->{locus}->{chr_start}-$crispr_well->{locus}->{chr_end}");
        push @output_data, \@output_row;
    }
}

my $csv = Text::CSV->new ({ binary => 1, auto_diag => 1 });
open my $fh, ">:encoding(utf8)", "pipeline1_data.csv" or die "pipeline1_data.csv: $!";
$csv->say ($fh, $_) for @output_data;
close $fh or die "pipeline1_data.csv: $!";


