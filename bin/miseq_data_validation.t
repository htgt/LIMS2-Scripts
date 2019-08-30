#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 79;

my @subs = qw(
                    miseq_data_validator
                    get_all_miseq_plates
                    get_resultset
                    check_plate
                    get_miseq_exps 
                    find_lost_wells 
                    get_well_ids 
                    get_well_name 
                    check_for_lost_well 
                    construct_file_path 
                    create_report_file 
                    write_report 
);

note('Testing package');

#BEGIN {
    use_ok( 'MiseqDataValidation', @subs );
#}

can_ok( __PACKAGE__, @subs );

note('Testing plate search');

my $schema = LIMS2::Model->new( user => 'lims2' )->schema;
my @test_plates = get_all_miseq_plates();
#ok ( scalar @test_plates > 1, 'More than 1 plate in resultset' );
foreach my $test_plate (@test_plates) {
    isa_ok( $test_plate, 'HASH' );
}
my $test_rs = get_resultset('Plate');
isa_ok( $test_rs, 'DBIx::Class::ResultSet');
my $test_miseq_10 = $test_rs->find({'name' => 'Miseq_010'})->as_hash;
my $test_miseq_11 = $test_rs->find({'name' => 'Miseq_011'})->as_hash;
check_plate($test_miseq_10);
ok( -e 'Miseq_010_lost_wells.txt', 'Report file created' );
check_plate($test_miseq_11);
ok( ! -e 'Miseq_011_lost_wells.txt', 'Report file not created' );



