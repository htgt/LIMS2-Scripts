#! /usr/bin/perl

use LIMS2::Model;
use feature qw/ say /;
use strict;
use warnings;
use Try::Tiny;

=head1 copy_plates -- Copy DNA plates in LIMS2

Syntax:
C<copy_plates.pl plate_synonym_list.txt>

This script invokes the LIMS2 create_plate_by_copy method to create copies of
DNA plates.

The input is a text file. Each line contains two plate names separated by whitespace.
The first plate name is the plate to be copied, the second plate name is the 
name of the copied plate.

The input file may be annotated with comments introduced by the '#' character.
Anything following the second plate name on each line will be ignored.

DJP-S 4 March 2013
=cut

say 'LIMS2 plate copy';

my $plates_file = $ARGV[0] || die 'Usage: copy_plates.pl plates file_name';

print 'Reading list of plates to copy from ' . $plates_file . "...\n";

open( my $plates_fh, '<', $plates_file )
    or die "Can't open file $plates_file: $! \n";

my @lines = <$plates_fh>;
close $plates_fh;

# Need to put each line in an array row because the input plate keys are not unique
# so we cannot put those straight into a hash.
my @plate_list;

my $plates_to_copy;
while (@lines) {
    my $line = shift @lines;
    chomp $line;
    next if $line =~ /^\s*#/;
    my ( $original_plate_name, $copy_plate_name ) = split( /\s+/, $line );
    push @plate_list, { $original_plate_name => $copy_plate_name };
    $plates_to_copy++;
}

say $plates_to_copy . ' plates to be copied.';
my $this_user = $ENV{'USER'} . '@sanger.ac.uk';
say 'Database user is: ' . $this_user;

my $model = LIMS2::Model->new( { user => 'webapp', audit_user => $this_user } );

my $n;
my $error_count;
foreach my $plate_line (@plate_list) {
    my ( $from_plate_name, $to_plate_name ) = each %{$plate_line};
    say "\n" . 'Copying plate: ' . $from_plate_name . ' to: ' . $to_plate_name;
    $model->txn_do(
        sub {
            try {
                $model->create_plate_by_copy(
                    {   from_plate_name => $from_plate_name,
                        to_plate_name   => $to_plate_name,
                        created_by      => $this_user,
                    }
                );
                $n++;
            }
            catch {
                say 'Error encountered while copying plate: ' . $_;
                $model->txn_rollback;
                $error_count++;
            }

        }
    );
}
say 'Copied ' . $n . ' of ' . $plates_to_copy . ' plates succesfully.'             if $n;
say 'Encountered errors in ' . $error_count . ' of ' . $plates_to_copy . ' plates' if $error_count;
say 'done';
