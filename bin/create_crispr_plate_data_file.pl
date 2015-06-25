#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Getopt::Long;
use Log::Log4perl ':easy';
use LIMS2::Model;
use Pod::Usage;
use Text::CSV;
use IO::Handle;
use Iterator::Simple qw( iterator );
use feature qw(say);

my $model = LIMS2::Model->new( user => 'lims2' );

my $log_level = $WARN;
GetOptions(
    'help'    => sub { pod2usage( -verbose    => 1 ) },
    'man'     => sub { pod2usage( -verbose    => 2 ) },
    'debug'   => sub { $log_level = $DEBUG },
    'verbose' => sub { $log_level = $INFO },
    'trace'   => sub { $log_level = $TRACE },
) or pod2usage(2);

my $crispr_file = shift @ARGV;
pod2usage(2) unless $crispr_file;

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );

open ( my $fh, '<', $crispr_file ) or die( "Can not open $crispr_file " . $! );
my $csv = Text::CSV->new();
$csv->column_names( @{ $csv->getline( $fh ) } );

my $io_output = IO::Handle->new_from_fd( \*STDOUT, 'w' );
my $output_csv = Text::CSV->new( { eol => "\n" } );
$output_csv->print( $io_output, [ 'well_name', 'crispr_id' ] );

my $well_name_iterator = well_names();

while ( my $data = $csv->getline_hr( $fh ) ) {
    my $crisprs_rs = $model->schema->resultset('Crispr')->search(
        {
            seq => $data->{seq}
        }
    );

    if ( $crisprs_rs->count == 1 ) {
        $output_csv->print( $io_output, [ $well_name_iterator->next, $crisprs_rs->first->id ] );
    }
    else {
        while ( my $crispr = $crisprs_rs->next ) {
            WARN( 'Crispr info: ' . $crispr->id );
        }
        WARN( 'Sequence: ' . $data->{seq} );
        LOGDIE( 'Seq matches multiple crisprs: ' . $crisprs_rs->count );
    }
}

sub well_names {
    my $letter = 'A';
    my $number = 01;

    iterator {
        if ( $number > 12 ) {
            return if $letter eq 'H';
            $number = 01;
            $letter++;
        }

        return $letter . sprintf( "%02d", $number++ );
    }
}

__END__

=head1 NAME

create_crispr_plate_data_file.pl - create CRISPR plate data file

=head1 SYNOPSIS

  create_crispr_plate_data_file.pl [crispr_data_file] [options]

      --help            Display a brief help message
      --man             Display the manual page
      --debug           Debug output
      --verbose         Verbose output
      --trace           Trace output

      Must specify crispr_data_file, a csv file with a seq column

=head1 DESCRIPTION

Given a csv file with crispr sequences produce a new CSV file with well names and crispr ids.
This file can be used to create a CRISPR plate in LIMS2.

=cut
