#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Data::UUID;
use IPC::Run qw( run );
use Smart::Comments;

my $dsn = 'dbi:Pg:host=htgt-db;port=5441;dbname=lims2_sp12';


if ( $dsn =~ /^dbi:Pg:host=(?<host>\S+);port=(?<port>\d+);dbname=(?<dbname>\S+)$/ ) {
    ### $+{host}
    ### $+{port}
    ### $+{dbname}
}
else {
    die( 'No' );

}
my $host = 'htgt-db'; 
my $port = 5441;
my $template_db_name = 'lims2_sp12';

my $db_user = 'lims2';
my $db_pass = 'team87';
my $temp_db_name;


#_create_temp_db();
#_drop_temp_db();

sub _create_temp_db {
    # createdb --host htgt-db --port 5441 --username lims2 --owner lims2 --template lims2_sp12 lims2_test_template

    $temp_db_name = Data::UUID->new->create_str();
    my @createdb_cmd = (
        'createdb',
        '--host', $host,
        '--port', $port,
        '--username', $db_user,
        '--owner', 'lims2',
        '--template', $template_db_name, 
        '--no-password',
        $temp_db_name
    );

    print "Created temp db $temp_db_name\n";
    local $ENV{PGPASSWORD} = $db_pass; 
    my $err = "";
    run( \@createdb_cmd,
        '<', \undef,
        '2>', \$err,
    ) or die( 'createdb command failed: ' . $err );

    return; 
}

sub _drop_temp_db {
    # dropdb --host htgt-db --port 5441 --username lims2 lims2_test_template

    my @dropdb_cmd = (
        'dropdb',
        '--host', $host,
        '--port', $port,
        '--username', $db_user,
        '--no-password',
        $temp_db_name
    );

    print "Dropping db $temp_db_name\n";
    local $ENV{PGPASSWORD} = $db_pass; 
    my $err = "";
    run( \@dropdb_cmd,
        '<', \undef,
        '2>', \$err,
    ) or die( 'dropdb command failed: ' . $err );

}
