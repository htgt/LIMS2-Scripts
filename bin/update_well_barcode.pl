#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Try::Tiny;
use Pod::Usage;


my $script_user = $ENV{USER}.'@sanger.ac.uk';
my $model = LIMS2::Model->new( user => 'webapp', audit_user => $script_user );

GetOptions(
    'help'         => sub { pod2usage( -verbose => 1 ) },
    'man'          => sub { pod2usage( -verbose => 2 ) },
    'barcode=s'    => \my $barcode,
    'new_state=s'  => \my $new_state,
    'comment=s'    => \my $comment,
);

die "You must provide a barcode" unless $barcode;
my $params = { barcode => $barcode };

die "You must provide a new barcode state" unless $new_state;
$params->{new_state} = $new_state;

$params->{user} = $script_user;

$params->{comment} = $comment if $comment;

$model->txn_do(
    sub {
        try{
            $model->update_well_barcode($params);
        }
        catch {
            print "barcode update failed: $_";
            $model->txn_rollback;
        };
    }
);

=head1 NAME

update_well_barcode.pl - change the state of a barcode in LIMS2

=head1 SYNOPSIS

  delete_well.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --barcode         Barcode to update
      --new_state       New barcode state, e.g. in_freezer
      --comment         Comment about the change to be stored in barcode event

=head1 DESCRIPTION

Changes the state of a barcode and creates a barcode event describing the change.

Note that this does NOT do any plate versioning. If your change affects the layout of a plate
then you need to use LIMS2::Model::Util::BarcodeActions.

=head1 AUTHOR

Anna Farne

=cut