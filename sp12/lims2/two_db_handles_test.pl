#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model::DBConnect;
use LIMS2::Model;
use Smart::Comments;

$ENV{LIMS2_DB} = 'LIMS2_STAGING'; 
my $schema_staging = LIMS2::Model::DBConnect->connect( 'LIMS2_DB', 'tasks' );
$ENV{LIMS2_DB} = 'LIMS2_TEST'; 
LIMS2::Model::DBConnect->_clear_connectors;
my $schema_test = LIMS2::Model::DBConnect->connect( 'LIMS2_DB', 'tasks' );

# pass in pre built schema object into model construction
my $model_staging = LIMS2::Model->new( user => 'tasks', schema => $schema_staging );
my $model_test = LIMS2::Model->new( user => 'tasks', schema => $schema_test );

test_model_db( $model_staging, 'LIMS2_STAGING' );
test_model_db( $model_test, 'LIMS2_TEST' );

sub test_model_db {
    my ( $model , $db ) = @_;

    my $dbh = $model->schema->storage->dbh;
    my $sth = $dbh->prepare_cached( 'SELECT current_catalog;' );
    $sth->execute( );
    my $data = $sth->fetchall_arrayref;
    ### $db
    ### $data
}
