#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Getopt::Long;
use Log::Log4perl ':easy';
use Pod::Usage;
use Text::CSV;
use YAML::Any;

my $log_level = $WARN;
GetOptions(
    'help'          => sub { pod2usage( -verbose => 1 ) },
    'man'           => sub { pod2usage( -verbose => 2 ) },
    'file=s'        => \my $file,
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );

LOGDIE( 'Must specify a input file' ) unless $file;

my $input_csv = Text::CSV->new();
open ( my $input_fh, '<', $file ) or die( "Can not open $file " . $! );
$input_csv->column_names( @{ $input_csv->getline( $input_fh ) } );
while ( my $data = $input_csv->getline_hr( $input_fh ) ) {
    print YAML::Any::Dump( $data );
}
close $input_fh;

__END__

=head1 NAME

csv_to_yaml.pl - convert a csv file to a yaml file 

=head1 SYNOPSIS

  csv_to_yaml.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --file            Input csv file

=head1 DESCRIPTION

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
