#!/usr/bin/env perl

# To run with DBIX trace on:
# DBIC_TRACE=1 perl -w -I ~/workspace/LIMS2-WebApp/lib ~/Sandbox/summaries_data_unload/summaries_test_data_unload.pl --plate_name XYZ00001 --well_name A01 --source_db LIMS2_LIVE --dest_db LIMS2_MYTESTDB

use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use LIMS2::Model::DBConnect;
use List::MoreUtils qw(uniq);
use Try::Tiny;                              # Exception handling
use Const::Fast;                            # Constant variables
use Getopt::Long;                           # Command line options
use Log::Log4perl ':easy';                  # DEBUG to INFO to WARN to ERROR to LOGDIE
use Smart::Comments;

#------------------------------------------------------------------
#  Variables
#------------------------------------------------------------------
my $origin_plate_name;                 # Origin plate name
my $origin_well_name;                  # Origin well name
my $loglevel = $INFO;                  # Logging level
my $schema_origin;                     # Schema for DB selecting from
my $schema_destination;                # Schema for DB inserting into
my $model_origin;                      # For connection to DB selecting from
my $model_destination;                 # For connection to DB inserting into
my $source_db;                         # source database accessor      e.g. LIMS2_LIVE
my $dest_db;                           # destination database accessor e.g. LIMS2_MYTESTDB

my %well_process_aux_data = (
    'create_di'              => \&_create_process_aux_data_create_di,
    'create_crispr'          => \&_create_process_aux_data_create_crispr,
    'int_recom'              => \&_create_process_aux_data_int_recom,
    '2w_gateway'             => \&_create_process_aux_data_2w_gateway,
    '3w_gateway'             => \&_create_process_aux_data_3w_gateway,
    'legacy_gateway'         => \&_create_process_aux_data_legacy_gateway,
    'final_pick'             => \&_create_process_aux_data_final_pick,
    'recombinase'            => \&_create_process_aux_data_recombinase,
    'cre_bac_recom'          => \&_create_process_aux_data_cre_bac_recom,
    'rearray'                => \&_create_process_aux_data_rearray,
    'dna_prep'               => \&_create_process_aux_data_dna_prep,
    'clone_pick'             => \&_create_process_aux_data_clone_pick,
    'clone_pool'             => \&_create_process_aux_data_clone_pool,
    'first_electroporation'  => \&_create_process_aux_data_first_electroporation,
    'second_electroporation' => \&_create_process_aux_data_second_electroporation,
    'freeze'                 => \&_create_process_aux_data_freeze,
    'xep_pool'               => \&_create_process_aux_data_xep_pool,
    'dist_qc'                => \&_create_process_aux_data_dist_qc,
);

GetOptions(
#    'help'             => sub { pod2usage( -verbose => 1 ) },
#    'man'              => sub { pod2usage( -verbose => 2 ) },
    'debug'            => sub { $loglevel = $DEBUG },
    'plate_name=s'     => \$origin_plate_name, # plate name to be copied across
    'well_name:s'      => \$origin_well_name,  # well name is optional
    'source_db'        => \$source_db,       # source database to copy from
    'dest_db'          => \$dest_db,         # destination database to copy into
);

# initialise logging
#Log::Log4perl->easy_init( { level => $loglevel, layout => '%p %m%n%x ' } );
Log::Log4perl->easy_init( { level => $loglevel, layout => '%p %m%n ' } );

# ---------------------------------
# Check input parameters
# ---------------------------------
if ( $origin_plate_name ) {
    INFO('Origin plate name : ' . $origin_plate_name);
}
else {
    ERROR('Origin plate name not supplied, cannot continue';
    die;
}

if ( $origin_well_name ) {
    INFO( 'Origin well name : ' . $origin_well_name );
}

if ( $source_db ) {
    INFO('Source DB : ' . $source_db);
}
else {
    ERROR('Source DB not supplied, cannot continue';
    die;
}

if ( $dest_db ) {
    INFO('Destination DB : ' . $dest_db);
}
else {
    ERROR('Destination DB not supplied, cannot continue';
    die;
}

# ---------------------------------
# Create new connections to DBs
# ---------------------------------
local $ENV{ LIMS2_DB } = $dest_db;
$schema_destination = LIMS2::Model::DBConnect->connect( 'LIMS2_DB', 'lims2' );
$model_destination = LIMS2::Model->new( user => 'lims2', schema => $schema_destination );

local $ENV{ LIMS2_DB } = $source_db;
# to allow setup of two database connections we need to clear connectors between connections
## no critic(ProtectPrivateSubs)
LIMS2::Model::DBConnect->_clear_connectors;
## use critic

$schema_origin = LIMS2::Model::DBConnect->connect( 'LIMS2_DB', 'lims2' );
$model_origin = LIMS2::Model->new( user => 'lims2', schema => $schema_origin );

# ---------------------------------
# Start copy process
# ---------------------------------
my ( $total_attempted, $total_unloaded ) = ( 0, 0 );

my $plate_orig = $model_origin->schema->resultset( 'Plate' )->find( { 'name' => $origin_plate_name } );

if ( ! $plate_orig ) {
    ERROR( "Failed to retrieve plate $origin_plate_name" );
    exit;
}

my ( $attempted, $unloaded ) = copy_plate_to_destination_db( $plate_orig, $origin_well_name );
$total_attempted += $attempted;
$total_unloaded  += $unloaded;

INFO( "Unloaded $total_unloaded of $total_attempted wells" );


sub copy_plate_to_destination_db {

    my ( $input_plate, $input_well_name ) = @_;

    #Log::Log4perl::NDC->push( $input_plate->name );

    my $plate_name = $input_plate->name;

    INFO 'Copying plate: ' . $plate_name . ' and well ' . $input_well_name;

    my ( $attempted, $unloaded ) = ( 0, 0 );

    my $plate_dest = try {
        retrieve_or_create_plate( $input_plate ) # in test db
    }
    catch {
        ERROR('Unable to retrieve or create plate ' . $_);
        undef;
    };

    #DEBUG 'Plate retrieved name : ' . $plate_dest->name if $plate_dest;

    return( $attempted, $unloaded ) unless $plate_dest;

    # loop through wells in the plate
    for my $curr_well ( $input_plate->wells ) {

        # match the well name input if supplied
        next if $input_well_name && $input_well_name ne $curr_well->name;

        # todo : anything special depending on type of design
        #if ($design->{type} eq 'conditional') {}
        $attempted++;

        #Log::Log4perl::NDC->push( $curr_well->name );
        try {
            copy_well_to_destination_db( $plate_dest, $curr_well );
            $unloaded++;
        }
        catch {
            ERROR( $_ );
        };
        #Log::Log4perl::NDC->pop;
    }
    #Log::Log4perl::NDC->pop;

    return ( $attempted, $unloaded );
}

sub copy_well_to_destination_db {

    my ( $plate_dest, $well_orig ) = @_;

    INFO "Copying well $well_orig";

    my $well_dest = retrieve_destination_well( $well_orig );

    if ( defined $well_dest ) {
        return $well_dest;
    }

    # Find the parent wells (if any, typically 1 e.g. XEP has many), returns ArrayRef
    my $parent_wells = find_parent_wells( $well_orig );

    if ( @{ $parent_wells } ) {
	    for my $parent_well ( @$parent_wells ) {
	    	# recursively set up parent wells first if they don't exist
	        copy_well_to_destination_db( retrieve_or_create_plate( $parent_well->plate ), $parent_well );
	    }
    }

    # fetch the process and plate types
    my $well_orig_process_type = $well_orig->output_processes->first->type_id;
    my $well_orig_plate_type   = $well_orig->plate->type_id;

    # construct well data for creation of well and processes
    my $well_data = build_well_data( $well_orig, $parent_wells, $well_orig_process_type );

    # If this is a create_di process, ensure the design exists before attempting to create the well
    if ( $well_orig_process_type eq 'create_di' ) {

        my $design_id = $well_data->{ 'process_data' }{ 'design_id' };

    	INFO 'Retrieving design for design ID : ' . $design_id;

        my $design = retrieve_or_create_design( $design_id );
    }

    INFO( "Copying well : $well_orig, process $well_orig_process_type, plate type $well_orig_plate_type" );

    # create the well
    $well_orig = create_well_in_destination_db( $well_data );

    return $well_orig;
}

sub retrieve_or_create_plate {
    my $plate = shift;
    return retrieve_destination_plate( $plate ) || create_destination_plate( build_plate_data($plate) );
}

sub retrieve_destination_plate {
    my $input_plate = shift;

    my $input_plate_name = $input_plate->name;

    DEBUG 'Retrieving plate ' . $input_plate_name;

    my $plate_dest = try {
    	$model_destination->schema->resultset( 'Plate' )->find( { 'name' => $input_plate_name } );
    }
    catch {
        $_->throw() unless $_->not_found;
        undef;
    };

    # check plate type
    if ( $plate_dest ) {
        my $plate_type = $plate_dest->type;

        DEBUG 'Plate ' . $input_plate_name . ' already exists in the destination DB, with plate type : ' . $plate_type if $plate_type;
    }

    return $plate_dest;
}

sub create_destination_plate {
    my $plate_data = shift;

    INFO( "Attempting to create plate $plate_data->{name}, type $plate_data->{type}" );
    ### $plate_data

    my $plate_dest = $model_destination->create_plate( $plate_data );

    return $plate_dest;
}

sub build_plate_data {
    my $plate = shift;

	# name        => { validate => 'plate_name' },
	# species     => { validate => 'existing_species', rename => 'species_id' },
	# type        => { validate => 'existing_plate_type', rename => 'type_id' },
	# description => { validate => 'non_empty_string', optional => 1 },
	# created_by => {
	# validate    => 'existing_user',
	# post_filter => 'user_id_for',
	# rename      => 'created_by_id'
	# },
	# created_at => { validate => 'date_time', optional => 1, post_filter => 'parse_date_time' },
	# comments   => { optional => 1 },
	# wells      => { optional => 1 },
	# is_virtual => { validate => 'boolean', optional => 1 },

    my %plate_data = (
        'name'        => $plate->name,
        'type'        => $plate->type_id,
        'description' => $plate->description,
        'created_by'  => 'test_user@example.org',
        'species'     =>  $plate->species_id,
    );

    return \%plate_data;
}


sub retrieve_destination_well {
    my $well = shift;

    my $plate_name = $well->plate->name;
    my $well_name = $well->name;

    my $well_dest = $model_destination->schema->resultset('Well')->find(
        {   'plate.name' => $plate_name,
            name         => $well_name
        },
        { join => 'plate', }
    );

    DEBUG "Well $well_dest already exists in the destination DB" if $well_dest;

    return $well_dest;
}

sub create_well_in_destination_db {
    my $well_data = shift;

    INFO( "Creating well $well_data->{plate_name}\[$well_data->{well_name}\], process type $well_data->{process_data}{type}" );

    my $well_dest = $model_destination->create_well( $well_data );

    return $well_dest;
}

sub build_well_data {
    my ( $input_well, $input_parent_wells, $input_well_process_type ) = @_;

    # plate_name   => { validate => 'existing_plate_name' },
    # well_name    => { validate => 'well_name', rename => 'name' },
    # accepted     => { validate => 'boolean', optional => 1 },
    # process_data => { validate => 'hashref' },
    # created_by => {
    #     validate    => 'existing_user',
    #     post_filter => 'user_id_for',
    #     rename      => 'created_by_id'
    # },
    # created_at => { validate => 'date_time', optional => 1, post_filter => 'parse_date_time' },

    my $plate = $input_well->plate;

    my %well_data = (
        'plate_name' => $plate->name,
        'well_name'  => $input_well->name,
        'created_by' => 'test_user@example.org',
#        'created_at' => $plate->created_at,
    );

    # fetch details from parent wells for process data
    my @parent_well_details;
    if ( $input_parent_wells ) {
    	for my $parent_well ( @$input_parent_wells ) {
    		push @parent_well_details, { 'plate_name' => $parent_well->plate->name, 'well_name' => $parent_well->name };
    	}
    }

    my %process_data = (
        'type'        => $input_well_process_type,
        #'input_wells' => @parent_well_details ? \@parent_well_details : []
    );

    if ( @parent_well_details ) {
        $process_data { 'input_wells' } = \@parent_well_details;
    }

    DEBUG "Setting up process data for well $input_well and process type $input_well_process_type";

    $well_process_aux_data{ $input_well_process_type }->( $input_well, \%process_data );

    DEBUG 'Process data created:' . (join ",", keys %process_data);

    $well_data{ 'process_data' } = \%process_data;

    return \%well_data;
}

sub _create_process_aux_data_create_di {
    my ( $well, $process_data ) = @_;

    DEBUG "create di Well $well";
    ### $process_data

    my $process = $well->output_processes->first;

    if ( $process->process_design ) {
        $process_data->{ 'design_id' } = $process->process_design->design->id;
    }

    #TODO put in bacs sp12 Fri 26 Jul 2013 13:28:27 BST
    #create_or_retrieve_bacs
    #$process_data { 'bacs' } = $well->bacs;

    return;
}

sub _create_process_aux_data_create_crispr{
    my ( $well, $process_data ) = @_;

    # no aux data

    return;
}

sub _create_process_aux_data_int_recom{
    my ( $well, $process_data ) = @_;

    DEBUG "Well $well";
    ### $process_data

    $process_data->{ 'cassette' } = $well->cassette->name;
    $process_data->{ 'backbone' } = $well->backbone->name;

    return;
}

sub _create_process_aux_data_2w_gateway{
    my ( $well, $process_data ) = @_;

    DEBUG "Well $well";
    ### $process_data

    #$process_data->{ 'cassette' } = $well->cassette->name;
    #$process_data->{ 'backbone' } = $well->backbone->name;

    my $process = $well->output_processes->first;

    
    if ( my $process_cassette = $process->process_cassette ) {
        $process_data->{ 'cassette' } = $process_cassette->name;
    }
    if ( my $process_backbone = $process->process_backbone ) {
        $process_data->{ 'backbone' } = $process_backbone->name;
    }

    my @recombinases = $process->process_recombinases->search( {} , { order_by => 'rank' } );

    if (@recombinases) {
        $process_data->{ 'recombinase' } = [ map { $_->recombinase_id } @recombinases ];
    }

    return;
}

sub _create_process_aux_data_3w_gateway{
    my ( $well, $process_data ) = @_;

    DEBUG "Well $well";
    ### $process_data

     $process_data->{ 'cassette' } = $well->cassette->name;
     $process_data->{ 'backbone' } = $well->backbone->name;

     my $process = $well->output_processes->first;

     my @recombinases = $process->process_recombinases->search( {} , { order_by => 'rank' } );

     if (@recombinases) {
         $process_data->{ 'recombinase' } = [ map { $_->recombinase_id } @recombinases ];
     }

    return;
}

sub _create_process_aux_data_legacy_gateway {
    my ( $well, $process_data ) = @_;

    DEBUG "Well $well";
    ### $process_data

    $process_data->{ 'cassette' } = $well->cassette->name;
    $process_data->{ 'backbone' } = $well->backbone->name;

    my $process = $well->output_processes->first;

    my @recombinases = $process->process_recombinases->search( {} , { order_by => 'rank' } );

    if (@recombinases) {
     $process_data->{ 'recombinase' } = [ map { $_->recombinase_id } @recombinases ];
    }

    return;
}

sub _create_process_aux_data_final_pick{
    my ( $well, $process_data ) = @_;

    #DEBUG "Well $well";
    ### $process_data

    # no aux data

    return;
}

sub _create_process_aux_data_recombinase{
    my ( $well, $process_data ) = @_;

    my $process = $well->output_processes->first;

    my @recombinases = $process->process_recombinases->search( {} , { order_by => 'rank' } );

    if (@recombinases) {
        $process_data->{ 'recombinase' } = [ map { $_->recombinase_id } @recombinases ];
    }

    return;
}

sub _create_process_aux_data_cre_bac_recom{
    my ( $well, $process_data ) = @_;

    # no aux data

    return;
}

sub _create_process_aux_data_rearray{
    my ( $well, $process_data ) = @_;

    # no aux data

    return;
}

sub _create_process_aux_data_dna_prep{
    my ( $well, $process_data ) = @_;

    # no aux data

    return;
}

sub _create_process_aux_data_clone_pick{
    my ( $well, $process_data ) = @_;

    # no aux data

    return;
}

sub _create_process_aux_data_clone_pool{
    my ( $well, $process_data ) = @_;

    # no aux data

    return;
}

sub _create_process_aux_data_first_electroporation{
    my ( $well, $process_data ) = @_;

    DEBUG "Attempting to get FEP cell line for $well";

    my $cell_line = $well->first_cell_line;

    DEBUG 'FEP cell line : ' . $cell_line->name if $cell_line;;

    die "No first_cell_line set for $well" unless $cell_line;

    $process_data->{ 'cell_line' } = $cell_line->name;

    return;
}

sub _create_process_aux_data_second_electroporation{
    my ( $well, $process_data ) = @_;

    # no aux data

    return;
}

sub _create_process_aux_data_freeze{
    my ( $well, $process_data ) = @_;

    # no aux data

    return;
}

sub _create_process_aux_data_xep_pool{
    my ( $well, $process_data ) = @_;

    # no aux data

    return;
}

sub _create_process_aux_data_dist_qc{
    my ( $well, $process_data ) = @_;

    # no aux data

    return;
}

sub retrieve_or_create_design {

    my $design_id = shift;

    DEBUG( 'In retrieve_or_create_design, design_id : ' . $design_id );

    return retrieve_destination_design( $design_id ) || create_destination_design( build_design_data( $design_id ) );
}

sub retrieve_destination_design {
    my $design_id = shift;

    INFO( "Attempting to retrieve design with ID :  $design_id" );

    my $design = try {
    	$model_destination->schema->resultset( 'Design' )->find( { 'id' => $design_id } );
    }
    catch {
    	$_->throw() unless $_->not_found;
        undef;
    };

    INFO( "Design $design->id already exists in the destination DB" ) if $design;

    return $design;
}

sub create_destination_design {
    my $design_data = shift;

    my $design_dest = $model_destination->create_design ( $design_data );

    INFO ('Design created has ID: ' . $design_dest->id ) if $design_dest;

    return $design_dest;
}

sub build_design_data {
    my $design_id = shift;

    DEBUG 'in build_design_data: design_id = ' . $design_id;

    Log::Log4perl::NDC->push( $design_id );

    my $design = $model_origin->schema->resultset( 'Design' )->find( { 'id' => $design_id } )
        or die "Failed to retrieve design with ID $design_id";

    my $design_data = get_design_data( $design );

    Log::Log4perl::NDC->pop;

    return $design_data;
}

sub get_design_data {
    my $design = shift;

    DEBUG( "in get_design_data" );

    # requirements as per Design plugin create method
    # my $design_data = {
    #     species                 => $design->species,
    #     id                      => $design->id,
    #     type                    => $design->type,
    #     created_at              => $design->created_at,
    #     created_by              => $design->created_by->id,
    #     phase                   => $design->phase,
    #     validated_by_annotation => $design->validated_by_annotation,
    #     name                    => $design->name,
    #     target_transcript       => $design->target_transcript,
    #     oligos                  => $design->oligos,
    #     genotyping_primers      => $design->genotyping_primers,
    #     comments                => $design->comments,
    #     gene_ids                => $design->gene_ids,
    # };

    # fetch data as hash from existing design object
    my $design_data = $design->as_hash( 1 );

    $design_data->{ gene_ids } = $design_data->{ assigned_genes };
    delete( $design_data->{ assigned_genes } );

    #each gene in gene_ids needs gene id and type id
    #gene_id      => $g->{gene_id},
    #gene_type_id => $g->{gene_type_id},

    $design_data->{ created_by } = 'test_user@example.org'; #$design->created_by->id;

    ### $design_data

    return $design_data;
}

sub find_parent_wells {

	my $child_well = shift;

    my @processes = $child_well->parent_processes;

    my $process = shift @processes;

    my @parent_wells = $process->input_wells->all;

    for my $well (@parent_wells) {
    	INFO "Parent well found : $well";
    }

    return \@parent_wells;
}
