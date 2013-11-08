#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Log::Log4perl ':easy';
use Pod::Usage;
use Try::Tiny;

my $log_level = $INFO;
GetOptions(
    'help'    => sub { pod2usage( -verbose    => 1 ) },
    'man'     => sub { pod2usage( -verbose    => 2 ) },
    'debug'   => sub { $log_level = $DEBUG },
    'plate=s' => \my $plate_name,
    'commit'  => \my $commit,
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %m%n' } );
my $model = LIMS2::Model->new( user => 'lims2' );

LOGDIE( "Must specify a plate with --plate option" ) unless $plate_name;

my $plate = $model->retrieve_plate( { name => $plate_name } );

LOGDIE( "$plate_name is not a design plate" ) unless $plate->type_id eq 'DESIGN';

$model->txn_do(
    sub {
        try{
            INFO( "Setting designs on $plate_name to cassette_first false" );
            for my $well( $plate->wells ) {
                my $design = $well->design;

                if ( $design->cassette_first ) {
                    DEBUG( "Design $design on well $well set to cassette_first false: " );
                    $design->update( { cassette_first => 0 } );
                }
                else {
                    WARN( "Design $design on well $well already cassette_first false: " );
                }
            }

            unless ( $commit ) {
                INFO( "Rollback" );
                $model->txn_rollback;
            }
        }
        catch {
            ERROR( "Error: $_" );
            $model->txn_rollback;
        };
    }
);

__END__

=head1 NAME

design_plate_set_cassette_first_false.pl - Set cassette_first flag as false on design plate

=head1 SYNOPSIS

  design_plate_set_cassette_first_false.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --debug           Debug output
      --plate           Name of design plate
      --commit          Commit changes, default is to rollback

=head1 DESCRIPTION

Get all the designs on a design plate and set the cassette_first flag to false.

=cut
