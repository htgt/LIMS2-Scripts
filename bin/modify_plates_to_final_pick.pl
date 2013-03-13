#! /usr/bin/perl

use LIMS2::Model;
use feature qw/ say /;
use strict;
use warnings;
use Try::Tiny;

=head1 modify_plates_to_final_pick

Syntax:
<modify_plates_to_final_pick.pl plate_list.txt>

This script runs SQL queries to update plates to FINAL_PICK

The input is a text file. Each line contains a plate name.
The plate name is the name of the plate to be converted to type FINAL_PICK.

The type of the process used to create each well in each plate will also be
updated depending on its current value:
'2w_' and '3w_gateway' will become 'legacy_gateway' 
'rearray' will become 'final_pick'

ACS 13 March 2013
=cut

say 'LIMS2 modify plates to final pick';

my $plates_file = $ARGV[0] || die 'Usage: modify_plates_to_final_pick.pl plates_list.txt';

print 'Reading list of plates to modify to final pick ' . $plates_file . "...\n";

open( my $plates_fh, '<', $plates_file )
    or die "Can't open file $plates_file: $! \n";

my @lines = <$plates_fh>;
close $plates_fh;

# Input plate names should be unique

# Connect to database via Model
my $this_user = $ENV{'USER'} . '@sanger.ac.uk';
say 'Database user is: ' . $this_user;

my $model = LIMS2::Model->new( { user => 'webapp', audit_user => $this_user } );

# counters
my $plates_processed = 0;
my $plates_updated = 0;
my $processes_checked = 0;
my $processes_updated_to_legacy = 0;
my $processes_updated_to_final_pick = 0;
my $error_count = 0;


$model->txn_do(
    sub{
        # cycle through the plates in the file
        while (@lines) {
            my $line = shift @lines;
            chomp $line;
            my $plate_name = $line;

            try {
                my $plate = $model->retrieve_plate({name=>$plate_name});
                $plates_processed++;
                print "Plate ID = ".$plate->id." and type = ".$plate->type_id."\n";

                # update the plate type to FINAL_PICK
                if($plate->type_id ne 'FINAL_PICK') {
                    print "Updating plate type to FINAL_PICK\n";
                    $plate->update({type_id=>'FINAL_PICK'});
					$plates_updated++;
                }
                # fetch the processes for the wells in the plate
                my @wells = $plate->wells;
                while(@wells) {
                    my $well = shift @wells;
                    print "Well ID = ".$well->id."\n";
                    my @processes= $well->parent_processes;
                    while(@processes) {
                        my $process = shift @processes;
                        $processes_checked++;
                        print "Process type = ".$process->type_id."\n";

                        if($process->type_id eq '2w_gateway' || $process->type_id eq '3w_gateway') {
                            print "Updating process type to legacy_gateway\n";
                            $process->update({type_id=>'legacy_gateway'});
                            $processes_updated_to_legacy++;
                        }
                        elsif($process->type_id eq 'rearray') {
                            print "Updating process type to final_pick\n";
                            $process->update({type_id=>'final_pick'});
                            $processes_updated_to_final_pick++;
                        }
                    }
                }
            } catch {
                say 'Error encountered while updating plate: ' . $plate_name;
                $error_count++;
            }
        }
    }
);
say 'Plates checked: '.$plates_processed;
say 'Plates updated: '.$plates_updated;
say 'Processes checked: '.$processes_checked;
say 'Processes updated to legacy_gateway: '.$processes_updated_to_legacy;
say 'Processes updated to final_pick: '.$processes_updated_to_final_pick;
say 'Encountered errors in '.$error_count.' plates' if $error_count;
say 'Done';
exit 0;