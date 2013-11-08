#! /usr/bin/perl

use LIMS2::Model;
use feature qw/ say /;
use LIMS2::Model::Util::CreateProcess qw( create_process_aux_data_recombinase );
use strict;
use warnings;

=head
For a given plate, which has a Dre recombinase on a clone_pick process,
remove the Dre recombinase from the process_recombinases.

Do this for all wells on the plate.

=cut

say "Remove Dre recombinase from plates";

my $this_user = $ENV{'USER'} . '@sanger.ac.uk';
say 'Database user is: ' . $this_user;

my $model = LIMS2::Model->new( { user => 'webapp', audit_user => $this_user } );

if ( !$ARGV[0] ) {
    say 'Usage: undre plate_name';
}

print "\n\nProcess types available:\n";
foreach my $proc_t ( @{ $model->list_process_types } ) {
    print $proc_t->id . ':' . $proc_t->description . "\n";
}
print "\n";

foreach my $plate_name ( @ARGV ) {
    modify_plate( $plate_name );
}

sub modify_plate {
    my $plate_name = shift;

    my $plate = $model->schema->resultset('Plate')->find( { name => $plate_name } );
    if ( !$plate ) {
        $model->throw( Validation => 'Plate ' . $plate_name . ' does not exist' . "\n" );
    }

    say 'Retrieved plate ' . $plate_name;

    my $well_rs = $model->schema->resultset('Well');

    my @wells_on_plate = $well_rs->search( { 'plate.name' => $plate_name }, { join => ['plate'], } );

    say 'Retrieved ' . @wells_on_plate . ' wells';

    foreach my $well (@wells_on_plate) {
        my @child_processes = $well->child_processes;
        say 'Well ' . $well->id . ' has ' . @child_processes . ' child processes';
        foreach my $child (@child_processes) {
            say '++ ' . $child->id . ' process_type: ' . $child->type_id;
            if ( $child->type_id eq 'clone_pick' ) {
                if ( $child->process_recombinases->first ) {
                    say '.. clone_pick: has recombinase ' . $child->process_recombinases->first->recombinase_id;

                    # update process type to rearray and delete the recombinase at this point
                    if ( $child->process_recombinases->first->recombinase_id eq 'Dre' ) {
                        remove_dre( $model, $child );
                    }
                }
            }
        }

    }
}


sub remove_dre {
    my $model = shift;
    my $child = shift;

    print "Dre will be removed\n";
    $model->txn_do(
        sub {
            $child->process_recombinases->first->delete;
        }
    );
    return;
}
