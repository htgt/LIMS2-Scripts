#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Try::Tiny;
use Pod::Usage;
use Log::Log4perl ':easy';

Log::Log4perl->easy_init( { level => $DEBUG } );

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

GetOptions(
    'help'         => sub { pod2usage( -verbose => 1 ) },
    'man'          => sub { pod2usage( -verbose => 2 ) },
    'plate_name=s' => \my $plate,
    'well_name=s'  => \my $well,
    'all!'         => \my $all,
    'commit'       => \my $commit,
);


# Find all assembly wells
my @assembly_wells;
if($plate){
	if($all){
        WARN "Plate name specified - ignoring 'all' flag";
	}
	if($well){
		my $well = $model->retrieve_well({
            plate_name => $plate,
            well_name => $well,
		});
		@assembly_wells = ($well);
	}
	else{
		my $plate = $model->retrieve_plate({
            name => $plate,
		});
		@assembly_wells = $plate->wells;
	}
}
elsif($all){
	my @assembly_plates = $model->schema->resultset('Plate')->search({
		    type_id => 'ASSEMBLY',
		    species_id => 'Human',
		})->all;
	@assembly_wells = map { $_->wells } @assembly_plates;
}

# Create hash with unique experiment key like design_id:nnnn-crispr_group_id:nnnnn
my $experiments = {};
foreach my $well (@assembly_wells){
    DEBUG "Generating experiment info for well $well";

    my $crispr_entity;
    try{
        $crispr_entity = $well->crispr_entity;
    };

    unless($crispr_entity){
    	WARN "No crispr entity found for well. Skipping experiment creation";
    	next;
    }

    my $design = $well->design;
    my ($gene_id) = $design->gene_ids;
    my $key = "design:".$design->id."-"
              .$crispr_entity->id_column_name.":".$crispr_entity->id;

    next if exists $experiments->{$key};

    my $experiment_info = {
    	design_id   => $design->id,
    	$crispr_entity->id_column_name => $crispr_entity->id,
    };

    # select project (species_id, gene_id). targeting_type = single_targeted
    my $project;
    try{
    	$project = $model->retrieve_project({
	        gene_id => $gene_id,
	        species_id => $well->plate->species_id,
	        targeting_type => 'single_targeted',
	    });
    };

    unless($project){
    	WARN "No project found for well $well. Skipping experiment creation";
    	next;
    }

    my $existing_experiment = $model->schema->resultset('Experiment')->search($experiment_info)->first;
    if($existing_experiment){
    	DEBUG "Experiment $key already exists in database";
    	$experiment_info->{exists} = 1;
    }

    DEBUG "Adding experiment $key to project ".$project->id;
    $experiment_info->{project_id} = $project->id;

    $experiments->{$key} = $experiment_info;
}

# Then create experiments
$model->txn_do(
    sub {
        try{
            foreach my $key (keys %$experiments){
            	my $experiment_info = $experiments->{$key};
            	next if $experiment_info->{exists};

                DEBUG "creating experiment $key";

                my $expt = $model->create_experiment($experiment_info);

                DEBUG "experiment created with ID ".$expt->id;
            }
            unless ( $commit ) {
                WARN "non-commit mode, rollback";
                $model->txn_rollback;
            }
        }
        catch {
            ERROR "create experiments failed: $_";
            $model->txn_rollback;
        };
    }
);

=head1 NAME

create_experiments_from_assembly_wells.pl - create entries in experiments table
for all assembly wells, or specified assembly well

=head1 SYNOPSIS

  delete_well.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --plate_name      Name of assembly plate
      --well_name       Name of assembly well
      --all             Or do it for all human assembly wells
      --commit          Commit changes, default is to rollback

=head1 DESCRIPTION

Create entries in the experiment table for all Human assembly wells,
or specified assembly well.

Assumes all assembly wells are for single_targeted projects.

=head1 AUTHOR

Anna Farne

=head1 BUGS

None reported... yet.

=head1 TODO

=cut