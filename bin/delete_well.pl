#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Try::Tiny;
use Pod::Usage;
use Perl6::Slurp;

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

GetOptions(
    'help'         => sub { pod2usage( -verbose => 1 ) },
    'man'          => sub { pod2usage( -verbose => 2 ) },
    'plate_name=s' => \my $plate,
    'well_name=s'  => \my $well,
    'barcode=s'    => \my $barcode,
    'commit'       => \my $commit,
);

if ( (!$well && !$plate) && !$barcode ) {
    die('Must specify plate and well, or barcode');
}

$model->txn_do(
    sub {
        try{
            my $params = {};
            my $name;
            if($barcode){
                $params = { barcode => $barcode };
                $name = "barcode: $barcode";
            }
            else{
                $params = { plate_name => $plate, well_name => $well };
                $name = "$plate [$well]";
            }
            $model->delete_well( $params );
            print "deleted well $name\n";
            unless ( $commit ) {
                print "non-commit mode, rollback\n";
                $model->txn_rollback;
            }
        }
        catch {
            print "delete well failed: $_\n";
            $model->txn_rollback;
        };
    }
);

=head1 NAME

delete_well.pl - delete a LIMS2 well

=head1 SYNOPSIS

  delete_well.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --plate_name      Name of plate well belongs to
      --well_name       Name of well to delete
      --commit          Commit changes, default is to rollback

=head1 DESCRIPTION

Delete a well in LIMS2, must specify plate the well belongs to.
Well deletion will not work if well has child wells.

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
