#!/usr/bin/env perl

##
# Selects PIQ clone data from LIMS2 into a csv file for use in updating iMits
#
# Used as a pre-cursor script in a jobrunner group, calling this first and sending the 
# resulting file into the htgt-batch/bin/update_piqs_in_tarmits.pl
##

# to run:
# perl -w ~/workspace/LIMS2-Scripts/bin/select_piq_clone_data.pl --lims2_data_filepath '/nfs/users/nfs_a/as28/tmp/piqs_to_imits/lims2_piq_data.csv'
# optionally add --debug  or add --mgi_gene_id MGI:123456

use strict;
use warnings FATAL => 'all';

use LIMS2::Model;                       # for LIMS2 database querying
use Log::Log4perl qw( :easy );          # for logging
use Getopt::Long;                       # for script input options
use Try::Tiny;                          # exception handling
use File::Basename qw( fileparse );     # split the entered filepath into directory and filename
use File::Path qw( make_path );         # to create a directory when it doesn't exist
use File::Spec;                         # to compile a filepath from directory and filename

my $mgi_gene_id;
my $lims2_data_filepath;
my $loglevel                    = $INFO;
my $LIMS2_DATA_FILENAME_DEFAULT = 'lims2_piq_data.csv';
my $counter_rows_written        = 0;

GetOptions(
    'mgi_gene_id:s'          => \$mgi_gene_id,
    'lims2_data_filepath:s'  => \$lims2_data_filepath,
    'debug'                  => sub { $loglevel = $DEBUG },
);

Log::Log4perl->easy_init(
	{
    	level  => $loglevel,
    	layout => '%p %m %n'
	}
);

if ( defined $mgi_gene_id ) { INFO 'Input gene MGI id : ' . $mgi_gene_id; }
if ( defined $lims2_data_filepath ) { INFO 'Input filepath <' . $lims2_data_filepath . '>' };

# check filepath input and die if not valid
unless ( defined $lims2_data_filepath ) {
    die 'ERROR: LIMS2 filepath not input to script, cannot continue';
}

# split the filepath
my ( $lims2_piq_data_file, $lims2_piq_data_directory ) = fileparse( $lims2_data_filepath );

DEBUG 'Filename = ' . $lims2_piq_data_file;
DEBUG 'File directory ' . $lims2_piq_data_directory;

# if the filepath does not include the filename assign the default
if ( !$lims2_piq_data_file ) {
    $lims2_piq_data_file = $LIMS2_DATA_FILENAME_DEFAULT;
    $lims2_data_filepath = File::Spec->catfile( $lims2_piq_data_directory, $lims2_piq_data_file );
}

# check the directory exists and create it if not
if ( !-d $lims2_piq_data_directory ) {
    make_path $lims2_piq_data_directory or die 'ERROR: Failed to create LIMS2 PIQ data directory: <' . $lims2_piq_data_directory . '>';
}

# backup the existing file if it exists
if ( -e $lims2_data_filepath ) {
    # file already exists
    DEBUG 'LIMS2 filepath <' . $lims2_data_filepath . '> already exists, creating backup';

    # backup existing file
    rename $lims2_data_filepath, $lims2_data_filepath . '.bak' or
                die 'ERROR: Could not rename existing lims2 PIQ data file <' . $lims2_data_filepath . "> and cannot continue. $!";
}

my @lims2_piq_data_rows;
_select_piq_data_from_lims2();

# write SQL to csv file at lims2_data_filepath, counting num rows
_write_csv_file( \@lims2_piq_data_rows );

# INFO message to state number of lines and filename written to
INFO 'LIMS2 data written to filepath <' . $lims2_data_filepath . '> successfully. Total rows: ' . $counter_rows_written;

##
# Select PIQ data from LIMS2
##
sub _select_piq_data_from_lims2 {
    # Create a new connection Model to link to DB
    my $model = LIMS2::Model->new( user => 'tasks' );

    # create SQL query and fetch results
    my $sql_query;
    if ( defined $mgi_gene_id ) {
        $sql_query = _create_sql_select_piqs_data_for_gene_from_lims2();
    }
    else {
        $sql_query = _create_sql_select_all_piqs_data_from_lims2();
    }
    my $sql_result = $model->schema->storage->dbh_do(
      sub {
         my ( $storage, $dbh ) = @_;
         my $sth = $dbh->prepare( $sql_query );
         $sth->execute or die "ERROR: Unable to execute query: $dbh->errstr\n";
         $sth->fetchall_arrayref({
         });
      }
    );

    # check SQL result valid
    unless ( scalar @$sql_result > 0 ) {
        ERROR 'Failed to select lims2 data from database';
        ERROR 'SQL:' . $sql_query;
        die 'ERROR: Failed to select LIMS2 PIQ data from database, cannot continue';
    }

    # parse SQL results into an array
    foreach my $row ( @$sql_result ) {
        my $curr_string = join ',', $row->{ 'project_name' }, $row->{ 'gene_mgi_id' }, $row->{ 'gene_marker_symbol' }, $row->{ 'clone_id' }, $row->{ 'piq_id' }, $row->{ 'piq_qc_pass' };
        push @lims2_piq_data_rows, $curr_string;
    }

    return;
}

##
# Write out to CSV logfile
##
sub _write_csv_file {
    my $lims2_piq_data_rows   = shift;

    unless ( scalar @$lims2_piq_data_rows > 0 ) { die 'No data in array sent to write to file, cannot continue'; }

    my $is_write_csv_ok = 1;

    try {
        my $csv_file = $lims2_data_filepath;

        if (-e $csv_file) {
            rename $csv_file, $csv_file . '.bak' or
                die 'ERROR: Could not rename lims2 data file ' . $csv_file . " $!";
        }

        my $OUTFILE;
        open $OUTFILE, '>>', $csv_file or
            die 'ERROR: Cannot open output lims2 data file ' . $csv_file . " for append: $!";

        foreach my $row ( @$lims2_piq_data_rows ) {
            print $OUTFILE ( $row )  . "\n";
            $counter_rows_written += 1;
        }
        close $OUTFILE;
    }
    catch {
        ERROR 'Failed to create lims2 data file at filepath '. $lims2_data_filepath;
        ERROR "Exception: ".$_;
        $is_write_csv_ok = 0;
    };

    unless ( $is_write_csv_ok ) { die 'ERROR: Could not write out to LIMS2 PIQ data file, cannot continue'; }

    return;
}

sub _create_sql_select_piqs_data_for_gene_from_lims2 {

    my $sql_query =  <<"SQL_END";
WITH cre_project_requests AS (
SELECT p.id AS project_id,
 p.htgt_project_id,
 p.sponsor_id,
 p.gene_id,
 p.targeting_type,
 pa.allele_type,
 pa.cassette_function,
 pa.mutation_type,
 cf.id AS cassette_function_id,
 cf.promoter,
 cf.conditional,
 cf.cre,
 cf.well_has_cre,
 cf.well_has_no_recombinase
FROM projects p
INNER JOIN project_alleles pa ON pa.project_id = p.id
INNER JOIN cassette_function cf ON cf.id = pa.cassette_function
WHERE p.sponsor_id   = 'Cre Knockin'
AND p.targeting_type = 'single_targeted'
AND p.species_id     = 'Mouse'
)
SELECT 'EUCOMMTOOLSCRE'                                AS project_name
, s.design_gene_id                                     AS gene_mgi_id
, s.design_gene_symbol                                 AS gene_marker_symbol
, (s.ep_pick_plate_name || '_' || s.ep_pick_well_name) AS clone_id
, (s.piq_plate_name || '_' || s.piq_well_name)         AS piq_id
, s.piq_well_accepted                                  AS piq_qc_pass
FROM summaries s
INNER JOIN cre_project_requests pr ON s.design_gene_id = pr.gene_id
WHERE s.design_type IN (SELECT design_type FROM mutation_design_types WHERE mutation_id = pr.mutation_type)
AND (
    (pr.conditional IS NULL)
    OR
    (pr.conditional IS NOT NULL AND s.final_pick_cassette_conditional = pr.conditional)
)
AND (
    (pr.promoter IS NULL)
    OR
    (pr.promoter IS NOT NULL AND pr.promoter = s.final_pick_cassette_promoter)
)
AND (
    (pr.cre IS NULL)
    OR
    (pr.cre IS NOT NULL AND s.final_pick_cassette_cre = pr.cre)
)
AND (
    (pr.well_has_cre IS NULL)
    OR
    (
        (pr.well_has_cre = true AND s.final_pick_recombinase_id = 'Cre')
        OR
        (pr.well_has_cre = false AND (s.final_pick_recombinase_id = '' OR s.final_pick_recombinase_id IS NULL))
    )
)
AND (
    (pr.well_has_no_recombinase IS NULL)
    OR
    (
     pr.well_has_no_recombinase IS NOT NULL AND (
      (pr.well_has_no_recombinase = true AND (s.final_pick_recombinase_id = '' OR s.final_pick_recombinase_id IS NULL))
       OR
      (pr.well_has_no_recombinase = false AND s.final_pick_recombinase_id IS NOT NULL)
     )
    )
)
AND s.piq_well_id > 0
AND s.design_gene_id       = '$mgi_gene_id'
GROUP BY s.design_gene_id
, s.design_gene_symbol
, (s.ep_pick_plate_name || '_' || s.ep_pick_well_name)
, (s.piq_plate_name || '_' || s.piq_well_name)
, s.piq_well_accepted
ORDER BY s.design_gene_symbol
, (s.ep_pick_plate_name || '_' || s.ep_pick_well_name)
, (s.piq_plate_name || '_' || s.piq_well_name)
SQL_END

    return $sql_query;

}

sub _create_sql_select_all_piqs_data_from_lims2 {

	my $sql_query =  <<"SQL_END";
WITH cre_project_requests AS (
SELECT p.id AS project_id,
 p.htgt_project_id,
 p.sponsor_id,
 p.gene_id,
 p.targeting_type,
 pa.allele_type,
 pa.cassette_function,
 pa.mutation_type,
 cf.id AS cassette_function_id,
 cf.promoter,
 cf.conditional,
 cf.cre,
 cf.well_has_cre,
 cf.well_has_no_recombinase
FROM projects p
INNER JOIN project_alleles pa ON pa.project_id = p.id
INNER JOIN cassette_function cf ON cf.id = pa.cassette_function
WHERE p.sponsor_id   = 'Cre Knockin'
AND p.targeting_type = 'single_targeted'
AND p.species_id     = 'Mouse'
)
SELECT 'EUCOMMTOOLSCRE'                                AS project_name
, s.design_gene_id                                     AS gene_mgi_id
, s.design_gene_symbol                                 AS gene_marker_symbol
, (s.ep_pick_plate_name || '_' || s.ep_pick_well_name) AS clone_id
, (s.piq_plate_name || '_' || s.piq_well_name)         AS piq_id
, s.piq_well_accepted                                  AS piq_qc_pass
FROM summaries s
INNER JOIN cre_project_requests pr ON s.design_gene_id = pr.gene_id
WHERE s.design_type IN (SELECT design_type FROM mutation_design_types WHERE mutation_id = pr.mutation_type)
AND (
    (pr.conditional IS NULL)
    OR
    (pr.conditional IS NOT NULL AND s.final_pick_cassette_conditional = pr.conditional)
)
AND (
    (pr.promoter IS NULL)
    OR
    (pr.promoter IS NOT NULL AND pr.promoter = s.final_pick_cassette_promoter)
)
AND (
    (pr.cre IS NULL)
    OR
    (pr.cre IS NOT NULL AND s.final_pick_cassette_cre = pr.cre)
)
AND (
    (pr.well_has_cre IS NULL)
    OR
    (
        (pr.well_has_cre = true AND s.final_pick_recombinase_id = 'Cre')
        OR
        (pr.well_has_cre = false AND (s.final_pick_recombinase_id = '' OR s.final_pick_recombinase_id IS NULL))
    )
)
AND (
    (pr.well_has_no_recombinase IS NULL)
    OR
    (
     pr.well_has_no_recombinase IS NOT NULL AND (
      (pr.well_has_no_recombinase = true AND (s.final_pick_recombinase_id = '' OR s.final_pick_recombinase_id IS NULL))
       OR
      (pr.well_has_no_recombinase = false AND s.final_pick_recombinase_id IS NOT NULL)
     )
    )
)
AND s.piq_well_id > 0
GROUP BY s.design_gene_id
, s.design_gene_symbol
, (s.ep_pick_plate_name || '_' || s.ep_pick_well_name)
, (s.piq_plate_name || '_' || s.piq_well_name)
, s.piq_well_accepted
ORDER BY s.design_gene_symbol
, (s.ep_pick_plate_name || '_' || s.ep_pick_well_name)
, (s.piq_plate_name || '_' || s.piq_well_name)
SQL_END
    return $sql_query;
}