#!/usr/bin/env perl

use strict;
use warnings;

use feature qw( say );
use Data::Dumper;

use LIMS2::Model;

use Text::CSV;
use autodie;
use Path::Class;

use Getopt::Long;
use Pod::Usage;

my ( @exon_ids, $species, $assembly, $csv_file );
GetOptions(
    'help'            => sub { pod2usage( 1 ) },
    'man'             => sub { pod2usage( 2 ) },
    'data-file=s'     => sub { my ( $name, $val ) = @_; $csv_file = file( $val ); },
    'species=s'       => sub { my ( $name, $val ) = @_; $species = ucfirst(lc $val); },
    'assembly=s'      => \$assembly,
) or pod2usage( 2 );

die "Please provide a species and an assembly" unless $species and $assembly;
die "Please provide a csv file with --data-file" unless $csv_file;

my $model = LIMS2::Model->new( user => 'lims2' );

my $fh = $csv_file->openr;

my $csv = Text::CSV->new();
$csv->column_names( @{ $csv->getline( $fh ) } );

my ( @pair_ids, $counter );
while ( my $data = $csv->getline_hr( $fh ) ) {
    my $id = $data->{exon_left_crispr_id} . "_" . $data->{exon_right_crispr_id};
    next if $id eq '_';

    push @pair_ids, $id;
    ++$counter;
}

say "Importing $counter pairs";

$model->import_wge_pairs( \@pair_ids, $species, $assembly );

1;

__END__

=head1 NAME

import_from_wge.pl -

=head1 SYNOPSIS

import_from_wge.pl [options]

    --help            Show this dialog
    --species         The species of the CRISPRs
    --assembly        The assembly of the CRISPRs
    --data-file       A CSV file containing WGE pair IDs

Example usage:

./import_from_wge.pl

=head1 DESCRIPTION



=head AUTHOR

Alex Hodgkins

=cut
