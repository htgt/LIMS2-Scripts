#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Pod::Usage;
use Try::Tiny;
use Perl6::Slurp;

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

GetOptions(
    'help'         => sub { pod2usage( -verbose => 1 ) },
    'man'          => sub { pod2usage( -verbose => 2 ) },
    'plate_name=s' => \my $plate,
    'commit'       => \my $commit,
);

my @plates;
if ( $plate ) {
    push @plates, $plate;
}
elsif ( $ARGV[0] ) {
    push @plates, map{ chomp; $_ } slurp( $ARGV[0] );
}
else {
    die( 'Must specify plate name(s)' );
}

for my $plate_name ( @plates ) {
    $model->txn_do(
        sub {
            try{
                $model->delete_plate( { name => $plate_name } );
                print "deleted plate $plate_name\n";
                unless ( $commit ) {
                    print "non-commit mode, rollback\n";
                    $model->txn_rollback;
                }
            }
            catch {
                print "delete plate failed: $_\n";
                $model->txn_rollback;
            };
        }
    );
}

=head1 NAME

delete_plate.pl - delete a LIMS2 plate

=head1 SYNOPSIS

  delete_plate.pl [options] input_file.txt

      --help            Display a brief help message
      --man             Display the manual page
      --plate_name      Name of plate to delete
      --commit          Commit changes, default is to rollback

=head1 DESCRIPTION

Delete a plate in LIMS2, either specify a plate name on the command line
or send in a file listing plate names seperated by newlines.

Deletion of a plate with child wells will not work.

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
