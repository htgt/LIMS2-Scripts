#!/usr/bin/env perl

use strict;
use warnings;

use LIMS2::Model;
use Try::Tiny;
use Getopt::Long;
use Log::Log4perl qw( :easy );
use Data::Dumper;

BEGIN { Log::Log4perl->easy_init( { level => $DEBUG, layout => '%p %m%n' } ) }

my %PROCESS_TYPE_DATA = (
    'int_recom'              => [qw( process_cassette process_backbone )],
    'cre_bac_recom'          => [qw( process_cassette process_backbone )],
    '2w_gateway'             => [qw( process_cassette process_backbone process_recombinases )],
    '3w_gateway'             => [qw( process_cassette process_backbone process_recombinases )],
    'recombinase'            => [qw( process_recombinases )],
    'clone_pick'             => [qw( process_recombinases )],
    'first_electroporation'  => [qw( process_cell_line process_recombinases )],
    'second_electroporation' => [qw( process_recombinases )],
    'crispr_ep'              => [qw( process_cell_line process_nuclease )],
    'crispr_vector'          => [qw( process_backbone )],
);

my %FIELD_NAMES = (
    process_backbone     => { relationship => "backbone",    column => "name" },
    process_cassette     => { relationship => "cassette",    column => "name" },
    process_recombinases => { relationship => "recombinase", column => "id"},
    process_cell_line    => { relationship => "cell_line",   column => "name"},
    process_nuclease     => { relationship => "nuclease",    column => "name" },
);

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

{
    my ( $target_well, $template_well ) = ( "A06", "A02" );
    my $user = "ah19\@sanger.ac.uk";

    # add_crispr_well( $user, "HCL0B", $target_well, 76878 );
    # add_crispr_well( $user, "HCR0B", $target_well, 76899 );

    # my @data = (
    #     #design
    #     {
    #         parent_plate  => "HG0",
    #         parent_well   => $target_well,
    #         target_plate  => "HINT0",
    #         target_well   => $target_well,
    #         template_well => $template_well,
    #     },
    #     {
    #         parent_plate  => "HINT0",
    #         parent_well   => $target_well,
    #         target_plate  => "HF0_Z",
    #         target_well   => $target_well,
    #         template_well => $template_well,
    #     },
    #     {
    #         parent_plate  => "HF0_Z",
    #         parent_well   => $target_well,
    #         target_plate  => "HFP0_Z",
    #         target_well   => $target_well,
    #         template_well => $template_well,
    #     },
    #     {
    #         parent_plate  => "HFP0_Z",
    #         parent_well   => $target_well,
    #         target_plate  => "HFPQ0_Z",
    #         target_well   => $target_well,
    #         template_well => $template_well,
    #     },
    #     #left_crispr
    #     {
    #         parent_plate  => "HCL0B",
    #         parent_well   => $target_well,
    #         target_plate  => "HCLS0B",
    #         target_well   => $target_well,
    #         template_well => $template_well,
    #     },
    #     {
    #         parent_plate  => "HCLS0B",
    #         parent_well   => $target_well,
    #         target_plate  => "HCLSQ0B",
    #         target_well   => $target_well,
    #         template_well => $template_well,
    #     },
    #     #right_crispr
    #     {
    #         parent_plate  => "HCR0B",
    #         parent_well   => $target_well,
    #         target_plate  => "HCRS0B",
    #         target_well   => $target_well,
    #         template_well => $template_well,
    #     },
    #     {
    #         parent_plate  => "HCRS0B",
    #         parent_well   => $target_well,
    #         target_plate  => "HCRSQ0B",
    #         target_well   => $target_well,
    #         template_well => $template_well,
    #     },
    #     #EP
    #     {
    #         parent_plate => "HFA0B",
    #         parent_well  => $target_well,
    #         target_plate => "HUEP0001B",
    #         target_well  => $target_well,
    #         template_well => $template_well,
    #     }
    # );

    for my $well ( @data ) {
        my $well = add_well_to_plate( 
            $well->{parent_plate},
            $well->{parent_well},
            $well->{target_plate},
            $well->{target_well},
            $well->{template_well},
            $user,
        );


    }

    #BRCA1 is in A08 on HFPQ0_Z

    # add_assembly_well( 
    #     'ah19@sanger.ac.uk', 
    #     "HFA0B", 
    #     $target_well, 
    #     [ 
    #         { plate_name => "HCLSQ0B", well_name => $target_well }, 
    #         { plate_name => "HCRSQ0B", well_name => $target_well }, 
    #         { plate_name => "HFPQ0_Z", well_name => "A07"} #design
    #     ] 
    # );
}

sub add_well_to_plate {
    my ( $parent_plate, $parent_well, $target_plate, $target_well, $template_well, $user ) = @_;

    my $db_template_well = $model->retrieve_well( {
        plate_name => $target_plate,
        well_name  => $template_well,
    } );

    DEBUG "Template well id is " . $db_template_well; 

    my $process = ($db_template_well->parent_processes)[0];

    DEBUG "Template process id is " . $process->id; 

    my %process_data = (
        type         => $process->type_id,
        input_wells  => [ { plate_name => $parent_plate, well_name => $parent_well } ],
        output_wells => [ { plate_name => $target_plate, well_name => $target_well } ] 
    );

    for my $field ( @{ $PROCESS_TYPE_DATA{$process->type_id} } ) {
        my @result = $process->$field;

        my $relationship = $FIELD_NAMES{$field}->{relationship};
        my $column       = $FIELD_NAMES{$field}->{column};

        #if more than one result then arrayref, otherwise single entry
        for my $entry ( @result ) {
            #if there's more than one entry for an item make it an array ref
            if ( defined $process_data{$relationship} ) {
                #item already exists, but is not an array, so make it one first.
                if ( ref $process_data{$relationship} ne 'ARRAY' ) {
                    $process_data{$relationship} = [ $process_data{$relationship} ];
                }

                #add new item to the array
                push @{ $process_data{$relationship} }, $entry->$relationship->$column;
            }
            else {
                $process_data{$relationship} = $entry->$relationship->$column;
            }
        }
    }

    print Dumper( \%process_data );

    return $model->txn_do(
        sub {
            my $well = $model->create_well( {
                plate_name   => $target_plate,
                well_name    => $target_well,
                process_data => \%process_data,
                created_by   => $user
            } );

            DEBUG "Created well " . $well;

            return $well;
        }
    );
}

#add a crispr to a top level crispr plate
#e.g. add_crispr_well( 'ah19@sanger.ac.uk', "HCR0", "A08", 74502 );
sub add_crispr_well {
    my ( $user, $plate_name, $well_name, $crispr_id ) = @_;

    return $model->txn_do(
        sub {
            return $model->create_well( {
                plate_name   => $plate_name, 
                well_name    => $well_name, 
                created_by   => $user,
                process_data => {
                    type         => "create_crispr", 
                    input_wells  => [], 
                    output_wells => [ { plate_name => $plate_name, well_name => $well_name} ], 
                    crispr_id    => $crispr_id,
                }, 
            } );
        }
    );
}

sub add_assembly_well {
    my ( $user, $plate_name, $well_name, $parent_plate_names ) = @_;

    return $model->txn_do(
        sub {
            return $model->create_well( {
                plate_name   => $plate_name, 
                well_name    => $well_name, 
                created_by   => $user,
                process_data => {
                    type         => "paired_crispr_assembly", 
                    input_wells  => $parent_plate_names, 
                    output_wells => [ { plate_name => $plate_name, well_name => $well_name } ]
                }, 
            } );
        }
    );
}

1;

=head1 NAME

add_well_to_plate.pl - delete a LIMS2 well

=head1 SYNOPSIS

  delete_well.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --plate_name      Name of plate well belongs to
      --well_name       Name of well to delete
      --commit          Commit changes, default is to rollback

=head1 DESCRIPTION

Delete a well in LIMS2, must specify plate the well belongs to.
Well deletion will not work if well has child wells.

=head1 AUTHOR

Sajith Perera
Alex Hodgkins

=cut