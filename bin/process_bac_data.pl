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
    'library=s'     => \my $bac_library,
    'assembly=s'    => \my $assembly,
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );

LOGDIE( 'Must specify a input file with --file' ) unless $file;
LOGDIE( 'Must specify a bac library with --library' ) unless $bac_library;
LOGDIE( 'Must specify a assembly with --assembly' ) unless $assembly;

my $input_csv = Text::CSV->new();
open ( my $input_fh, '<', $file ) or die( "Can not open $file " . $! );
$input_csv->column_names( @{ $input_csv->getline( $input_fh ) } );
while ( my $data = $input_csv->getline_hr( $input_fh ) ) {
    process_bac_data( $data );
}
close $input_fh;

sub process_bac_data {
    my ( $data ) = @_;

    print YAML::Any::Dump( 
        {
            bac_name => $data->{bac_name},
            bac_library => $bac_library,
            loci => [
                {
                    assembly => $assembly,
                    chr_name => $data->{chr_name},
                    chr_start => $data->{chr_start},
                    chr_end => $data->{chr_end},
                }
            ]
        }
    );
}

__END__

=head1 NAME

process_bac_data.pl - process bac data in csv file to a yaml file 

=head1 SYNOPSIS

  process_bac_data.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --library         Name of BAC library
      --assembly        Assembly of BAC locus information
      --file            Input csv file

=head1 DESCRIPTION

Take a csv file of bac clone data, with following field headers:
-bac_name
-chr_name
-chr_start
-chr_end

Output a yaml file that can be loaded into LIMS2 via load-bacs task.
YAML format

bac_name: RP11-23423
bac_library: RP11
loci:
   -assembly: GRCh37
    chr_name: X
    chr_start: 123423
    chr_end: 123787

=cut
