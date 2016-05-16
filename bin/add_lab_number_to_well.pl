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
    'well=s'       => \my $well_name,
    'plate=s'      => \my $plate_name,
    'lab_number=s' => \my $lab_number,
);

unless ($well_name and $plate_name and $lab_number){
	pod2usage( -verbose => 1);
	die;
}

my $well = $model->retrieve_well({ well_name => $well_name, plate_name => $plate_name });

die "Well $well_name $plate_name not found" unless $well;

$model->txn_do(
    sub {
        try{
            $model->create_well_lab_number({
                well_id    => $well->id,
                lab_number => $lab_number,
                created_by => $script_user,
            });
        }
        catch {
            print "Lab number creation failed: $_\n";
            $model->txn_rollback;
        };
    }
);

=head1 NAME

add_lab_number_to_well.pl - add a lab number to a well in LIMS2

=head1 SYNOPSIS

  add_lab_number_to_well.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --well            Name of well to update
      --plate           Name of plate to update
      --lab_number      Lab number to add

=head1 DESCRIPTION

Adds a lab number to an existing well in LIMS2 (lab numbers must be unique)

=head1 AUTHOR

Anna Farne

=cut
