#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Try::Tiny;
use Pod::Usage;
use feature qw( say );
use List::MoreUtils qw( uniq );

my $model = LIMS2::Model->new( user => 'lims2' );

GetOptions(
    'help'           => sub { pod2usage( -verbose => 1 ) },
    'man'            => sub { pod2usage( -verbose => 2 ) },
    'plate=s'        => \my $plate_name,
    'well=s'         => \my $well_name,
    'parent-plate=s' => \my $parent_plate_name,
    'parent-well=s'  => \my $parent_well_name,
    'commit'         => \my $commit,
);

my ( $well, $parent_well );
if ( $plate_name && $well_name ) {
    $well = $model->retrieve_well( { plate_name => $plate_name, well_name => $well_name } );
}
else {
    pod2usage('Must specify a plate and well name');
}

if ( $parent_plate_name && $parent_well_name ) {
    $parent_well = $model->retrieve_well(
        { plate_name => $parent_plate_name, well_name => $parent_well_name } );
}
else {
    pod2usage('Must specify a parent plate and parent well name');
}

say 'NOTE: This script keeps the same process and all its related data.';
say '      It only changes the input / parent well.';
say '      The only check done is to see if the current parent plate type is the same as the new one.';

say "Reparenting well $well";

# work out current parent well
my @process_input_wells = $well->process_output_wells->first->process->process_input_wells;
if ( scalar( @process_input_wells ) > 1 ) {
    die( "Well has multiple input wells, this script cannot handle this yet" );
}
my $process_input_well = shift @process_input_wells;
my $current_parent_well = $process_input_well->well;

say ".. current parent well is $current_parent_well";

my $parent_well_type = $parent_well->plate->type_id;
my $current_parent_well_type = $current_parent_well->plate->type_id;
if ( $parent_well_type ne $current_parent_well_type ) {
    die( "Current parent well type $current_parent_well_type does not match the new parent well type $parent_well_type" );
}

$model->txn_do(
    sub {
        try{
            say ".... new parent well will be $parent_well";
            $process_input_well->update( { well_id => $parent_well->id  } );
            unless ( $commit ) {
                print "non-commit mode, rollback\n";
                $model->txn_rollback;
            }
        }
        catch {
            print "reparent well failed: $_\n";
            $model->txn_rollback;
        };
    }
);

=head1 NAME

reparent_well.pl - Change the parent of a well, keeps the same process.

=head1 SYNOPSIS

  reparent_well.pl [options]
      --help            Display a brief help message
      --man             Display the manual page
      --plate           Name of plate
      --well            Name of well
      --parent-plate    Name of parent plate
      --parent-well     Name of parent well
      --commit          Commit the changes, default to rollback

The plate types of the current parent well and the new parent well must match.
This is because this script does not try to change the process, it keeps the same one.

=head1 DESCRIPTION

Use this script with caution, make sure you know exactly what you are doing.
Reparent a well in LIMS2, keeps the same process data, only the input well changes.
The script does a basic check to see the new and old parent wells are on plates of the same type.

Also only works for wells with a single parent well.

=cut
