#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Getopt::Long;
use Pod::Usage;
use IO::Handle;
use Text::CSV;
use YAML::Any qw( LoadFile );
use Perl6::Slurp;

GetOptions(
    'help' => sub { pod2usage( -verbose => 1 ) },
    'man'  => sub { pod2usage( -verbose => 2 ) },
) or pod2usage(2);

my $file = $ARGV[0];
die( 'Must specify a input file' ) unless $file;

my @data = LoadFile( $file );
my @column_names;
if ( my $col_file = $ARGV[1] ) {
    @column_names = slurp $col_file, {chomp=>1};
}
else {
    @column_names = keys %{ $data[0] };
}

my $io_output = IO::Handle->new_from_fd( \*STDOUT, 'w' );
my $output_csv = Text::CSV->new( { eol => "\n" } );
$output_csv->print( $io_output, \@column_names );

for my $record ( @data ) {
    my @array_data = @{ $record }{ @column_names };
    $output_csv->print( $io_output, \@array_data );
}

__END__

=head1 NAME

yaml_to_csv.pl - convert a yaml file ( array of hashes ) to a csv file

=head1 SYNOPSIS

  yaml_to_csv.pl [options] filename [column_filename]

      --help            Display a brief help message
      --man             Display the manual page

First argument should be a name of the yaml file.
Second optinal argument is a text file with the output csv column names, one to each line.

=head1 DESCRIPTION

Convert a YAML file, with data stored in array of hashes format only, to a csv.
The hash data will only be one level deep, anything more can not be handled.

=cut
