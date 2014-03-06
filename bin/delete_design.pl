#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Pod::Usage;
use Try::Tiny;

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

GetOptions(
    'help'     => sub { pod2usage( -verbose => 1 ) },
    'man'      => sub { pod2usage( -verbose => 2 ) },
    'design=i' => \my $design_id,
    'commit'   => \my $commit,
);

die( 'Must specify design id' ) unless $design_id;

$model->txn_do(
    sub {
        try{
            my $design = $model->c_retrieve_design( { id => $design_id } );
            $design->genes->delete;
            $model->c_delete_design( { id => $design_id, cascade => 1 } );
            print "deleted design $design_id\n";
            unless ( $commit ) {
                print "non-commit mode, rollback\n";
                $model->txn_rollback;
            }
        }
        catch {
            print "delete design failed: $_\n";
            $model->txn_rollback;
        };
    }
);

__END__

=head1 NAME

delete_design.pl - delete a LIMS2 design

=head1 SYNOPSIS

  delete_design.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --design          The id of the design
      --commit          Commit changes, default is to rollback

=head1 DESCRIPTION

Delete a LIMS2 design, will not work if design belongs to a design well.

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
