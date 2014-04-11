#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Pod::Usage;
use Perl6::Slurp;
use feature qw( say );
use List::MoreUtils qw( uniq );

my $model = LIMS2::Model->new( user => 'lims2' );

GetOptions(
    'help'      => sub { pod2usage( -verbose => 1 ) },
    'man'       => sub { pod2usage( -verbose => 2 ) },
    'plate=s'  => \my $plate_name,
);

my @plates;
if ( $plate_name ) {
    push @plates, $plate_name;
}
elsif ( $ARGV[0] ) {
    push @plates, map{ chomp; $_ } slurp( $ARGV[0] );
}
else {
    podusage(2);
}

my @design_ids;
for my $name ( @plates ) {
    my $plate = $model->retrieve_plate( { name => $name } );

    die "unable to locate plate $name" unless $plate;

    push @design_ids, map{ $_->design->id } $plate->wells->all;
}

say $_ for uniq @design_ids;

=head1 NAME

designs_from_plate.pl - Get all the designs linked to a plate

=head1 SYNOPSIS

  designs_from_plate.pl [options] [file]

      --help            Display a brief help message
      --man             Display the manual page
      --plate           Name of a plate

  You must specify either a single plate name through the command line
  or specify a file with a list of plate names, one per line.

=head1 DESCRIPTION

Output the designs linked to all the wells on a given plate or plates.

=cut
