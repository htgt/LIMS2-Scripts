#! /usr/bin/perl

use LIMS2::Model;
use feature qw/ say /;
use LIMS2::Model::Util::CreateProcess qw( create_process_aux_data_recombinase );
use strict;
use warnings;

=head
For a given plate, which has no dre recombinase on the first electroporation,
and where the plate name is CEP <= 26
get the process for the input well and make that a recombinase process with dre.

Do this for all wells on the plate.

=cut

say "Add Dre recombinase to certain plates";

my $this_user = $ENV{'USER'} . '@sanger.ac.uk';
say 'Database user is: ' . $this_user;

my $model = LIMS2::Model->new( { user => 'webapp', audit_user => $this_user } );

if ( !$ARGV[0] ) {
    say 'Usage: add_dre plate_name';
}

print "\n\nProcess types available:\n";
foreach my $proc_t ( @{ $model->list_process_types } ) {
    print $proc_t->id . ':' . $proc_t->description . "\n";
}
print "\n";

foreach my $plate_name ( @ARGV ) {
    modify_plate( $plate_name );
}

exit(1);

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

    foreach my $well ( @wells_on_plate ) {
       my @parent_processes = $well->parent_processes;
       say "\n" .'Well ' . $well->id . ' has ' . @parent_processes . ' parent processes';
       foreach my $parent ( @parent_processes ) {
            say '-- ' . $parent->id . ' process_type: ' . $parent->type_id;
            if ( $parent->process_recombinases->first ) {
                say '.. recombinase already exists: ' . $parent->process_recombinases->first->recombinase_id;
            }
            else {
                print ".. adding Dre recombinase\n";
                add_dre_to_FEP( $model, 'Dre', $parent );
            }
       }
    }
}

sub add_dre_to_FEP {
    my $model       = shift;
    my $recombinase = shift;
    my $fep_process = shift;

    if ( $fep_process->type_id eq 'first_electroporation' ) {
        $model->txn_do(
            sub {
                create_process_aux_data_recombinase( $model, { recombinase => [$recombinase] }, $fep_process );
            }
        );
    }

    return;
}
