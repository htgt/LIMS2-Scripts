#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Log::Log4perl ':easy';
use Pod::Usage;
use feature qw( say );
use Iterator::Simple qw( iterator );


my $log_level = $INFO;
GetOptions(
    'help'     => sub { pod2usage( -verbose    => 1 ) },
    'man'      => sub { pod2usage( -verbose    => 2 ) },
    'debug'    => sub { $log_level = $DEBUG },
    'plate=s'  => \my $plate_name,
    'force-96' => \my $force_96_well,
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %m%n' } );
my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

LOGDIE( 'Must specify a plate name' ) unless $plate_name;

my $plate = $model->retrieve_plate( { name => $plate_name } );

INFO( 'Plate: ' . $plate->name . ' of type ' . $plate->type_id);

say 'well_name,parent_plate,parent_well';
my @wells = $plate->wells->search(
    { }, { order_by => 'name' }
);

my $well_name_iterator = well_names();

MAIN: while ( my $well_name_96 = $well_name_iterator->next ) {
    my $well = shift @wells;
    while ( !$well || $well_name_96 ne $well->name ) {
        last MAIN unless $well_name_96;
        say "$well_name_96,-";
        $well_name_96 = $well_name_iterator->next;
    }
   my @process = $well->parent_processes;
   LOGDIE( 'More than one parent process??' ) if scalar( @process ) > 1;
   my $process = shift @process;

   my @input_wells = $process->input_wells;
   LOGDIE( 'More than one input wells??' ) if scalar( @input_wells ) > 1;
   my $input_well = shift @input_wells;

   INFO(
       'Well: ' . $well->name 
       . ' Input Plate: ' . $input_well->plate->name 
       . ' Input Well: ' . $input_well->name 
       . ' Input Plate Type: ' . $input_well->plate->type_id 
       . ' Process Type: ' . $process->type_id
   );
   say $well->name . ',' . $input_well->plate->name . ',' . $input_well->name;

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

plate_parent_wells.pl - Find parent well(s) for a given plates wells 

=head1 SYNOPSIS

  plate_parent_wells.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --debug           debug output

=head1 DESCRIPTION

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
