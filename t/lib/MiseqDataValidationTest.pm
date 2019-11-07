package MiseqDataValidationTest;

use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/../../lib";
use base qw( Test::Class );
use Test::More tests => 27;
use Test::Deep;
use Test::MockModule;
use LIMS2::Test model => { classname => __PACKAGE__ };
use MiseqDataValidation qw(:all);
use Data::Dumper;

sub redefined_construct_file_path {
    my ( $plate_name, $i, $miseq_exp_name ) = @_;
    my $file_path =
      "$Bin/../../test_fixtures/quant_files/$plate_name/$miseq_exp_name/$i.txt";
    return $file_path;
}

undef &construct_file_path;
*MiseqDataValidation::construct_file_path = \&redefined_construct_file_path;

sub test_module_import : Test(2) {
    my @subs = qw(
      miseq_data_validator
      get_all_miseq_plates
      check_plate
      get_mismatches
      get_miseq_exps
      find_lost_wells
      get_well_names
      check_for_lost_well
      construct_file_path
      create_report_file
      write_report
    );
    use_ok( 'MiseqDataValidation', @subs );
    can_ok( __PACKAGE__, @subs );

    return;
}

sub test_get_all_miseq_plates : Test(1) {
    my @test_plates     = get_all_miseq_plates();
    my @expected_plates = (
        {
            'id'          => 12307,
            'name'        => 'Miseq010',
            'description' => '',
            'type'        => 'MISEQ',
            'created_by'  => 'beth',
            'created_at'  => '2019-11-05T16:33:42',
            'wells'       => [ 'A02', 'C09', 'K09', 'O02' ]
        },
        {
            'id'          => 12316,
            'name'        => 'Miseq011',
            'description' => '',
            'type'        => 'MISEQ',
            'created_by'  => 'beth',
            'created_at'  => '2019-11-05T16:33:42',
            'wells'       => [
                '',    'A01', 'B04', 'C02', 'C14', 'C17',
                'E09', 'F04', 'I01', 'P09'
            ]
        },
        {
            'id'          => 13579,
            'name'        => 'Test_Plate',
            'description' => '',
            'type'        => 'MISEQ',
            'created_by'  => 'beth',
            'created_at'  => '2019-11-05T16:33:42',
            'wells'       => []
        }
    );
    cmp_deeply( \@test_plates, bag(@expected_plates),
        'Expected plate data retrieved' );

    return;
}

sub test_check_plate : Test(3) {
    my $test_miseq_10 = {
        'id'          => 12307,
        'name'        => 'Miseq010',
        'description' => '',
        'type'        => 'MISEQ',
        'created_by'  => 'beth',
        'created_at'  => '2019-11-05T16:33:42',
        'wells'       => [ 'A02', 'C09', 'K09', 'O02' ]
    };
    check_plate($test_miseq_10);
    ok( -e 'Miseq010_lost_wells.txt', 'Report file created for Miseq 10' );

    open( my $report_fh, '<', 'Miseq010_lost_wells.txt' )
      or die "Couldn't open Miseq 10 report file: $!";
    chomp( my @report_file_data = <$report_fh> );
    close $report_fh or die "Couldn't close Miseq 10 report file: $!";
    my @expected_file_data =
      ( 'CG_PDX_1_FPDM2: H05 M11', 'Plate_2: K01 O11', 'HUPEPD0017_1: M06' );
    cmp_deeply(
        \@report_file_data,
        bag(@expected_file_data),
        'Correct data in Miseq 10 report file'
    );

    my $test_miseq_11 = {
        'id'          => 12316,
        'name'        => 'Miseq011',
        'description' => '',
        'type'        => 'MISEQ',
        'created_by'  => 'beth',
        'created_at'  => '2019-11-05T16:33:42',
        'wells' =>
          [ '', 'A01', 'B04', 'C02', 'C14', 'C17', 'E09', 'F04', 'I01', 'P09' ]
    };
    check_plate($test_miseq_11);
    ok( !-e 'Miseq011_lost_wells.txt', 'Report file not created for Miseq 11' );

    return;
}

sub test_get_mismatches : Test(2) {
    my $miseq_10_mismatches          = get_mismatches( 12307, 'Miseq010' );
    my $expected_miseq_10_mismatches = {
        'CG_PDX_1_FPDM2' => [ 'H05', 'M11' ],
        'Plate_2'        => [ 'K01', 'O11' ],
        'HUPEPD0017_1'   => ['M06']
    };
    cmp_deeply( $miseq_10_mismatches, $expected_miseq_10_mismatches,
        'Correct mismatches for Miseq 10' );

    my $miseq_11_mismatches          = get_mismatches( 12316, 'Miseq011' );
    my $expected_miseq_11_mismatches = {};
    is_deeply( $miseq_11_mismatches, $expected_miseq_11_mismatches,
        'Correct mismatches for Miseq 11' );

    return;
}

sub test_get_miseq_exps : Test(4) {
    my @miseq_10_exps          = get_miseq_exps(12307);
    my @expected_miseq_10_exps = (
        {
            'id'              => 170,
            'miseq_id'        => 13,
            'name'            => 'CG_PDX_1_FPDM2',
            'experiment_id'   => undef,
            'parent_plate_id' => undef,
            'gene'            => undef,
            'nhej_count'      => undef,
            'read_count'      => undef
        },
        {
            'id'              => 165,
            'miseq_id'        => 13,
            'name'            => 'Plate_2',
            'experiment_id'   => undef,
            'parent_plate_id' => undef,
            'gene'            => undef,
            'nhej_count'      => undef,
            'read_count'      => undef
        },
        {
            'id'              => 167,
            'miseq_id'        => 13,
            'name'            => 'HUPEPD0017_1',
            'experiment_id'   => undef,
            'parent_plate_id' => undef,
            'gene'            => undef,
            'nhej_count'      => undef,
            'read_count'      => undef
        }
    );
    print Dumper( \@miseq_10_exps );
    cmp_deeply(
        \@miseq_10_exps,
        bag(@expected_miseq_10_exps),
        'Correct experiments for Miseq 10'
    );

    my @miseq_11_exps          = get_miseq_exps(12316);
    my @expected_miseq_11_exps = (
        {
            'id'              => 339,
            'miseq_id'        => 21,
            'name'            => 'DDD_SETD1B_2_CAC',
            'experiment_id'   => undef,
            'parent_plate_id' => undef,
            'gene'            => undef,
            'nhej_count'      => undef,
            'read_count'      => undef
        },
        {
            'id'              => 340,
            'miseq_id'        => 21,
            'name'            => 'HUEDQ0501_B_KMT2C_1',
            'experiment_id'   => undef,
            'parent_plate_id' => undef,
            'gene'            => undef,
            'nhej_count'      => undef,
            'read_count'      => undef
        },
        {
            'id'              => 288,
            'miseq_id'        => 21,
            'name'            => 'SETD5_2_PLATE11C',
            'experiment_id'   => undef,
            'parent_plate_id' => undef,
            'gene'            => undef,
            'nhej_count'      => undef,
            'read_count'      => undef
        }
    );
    cmp_deeply(
        \@miseq_11_exps,
        bag(@expected_miseq_11_exps),
        'Correct experiments for Miseq 11'
    );

    my @no_miseq_exps          = get_miseq_exps(13579);
    my @expected_no_miseq_exps = ();
    is_deeply( \@no_miseq_exps, \@expected_no_miseq_exps,
        'Returns empty list if no Miseq experiments linked to plate' );

    my @no_miseq_plate          = get_miseq_exps(54321);
    my @expected_no_miseq_plate = ();
    is_deeply( \@no_miseq_plate, \@expected_no_miseq_plate,
        'Returns empty list if no Miseq plate found' );

    return;
}

sub test_find_lost_wells : Test(3) {
    my $test_exp_1 = {
        'id'              => 165,
        'miseq_id'        => 13,
        'name'            => 'Plate_2',
        'experiment_id'   => undef,
        'parent_plate_id' => undef,
        'gene'            => undef,
        'nhej_count'      => undef,
        'read_count'      => undef
    };
    my @exp_1_lost_wells          = find_lost_wells( $test_exp_1, 'Miseq010' );
    my @expected_exp_1_lost_wells = ( 'K01', 'O11' );
    cmp_deeply(
        \@exp_1_lost_wells,
        bag(@expected_exp_1_lost_wells),
        'Correct lost wells returned for test experiment 1'
    );

    my $test_exp_2 = {
        'id'              => 167,
        'miseq_id'        => 13,
        'name'            => 'HUPEPD0017_1',
        'experiment_id'   => undef,
        'parent_plate_id' => undef,
        'gene'            => undef,
        'nhej_count'      => undef,
        'read_count'      => undef
    };
    my @exp_2_lost_wells          = find_lost_wells( $test_exp_2, 'Miseq010' );
    my @expected_exp_2_lost_wells = ('M06');
    is_deeply( \@exp_2_lost_wells, \@expected_exp_2_lost_wells,
        'Correct lost well returned for test experiment 2' );

    my $test_exp_3 = {
        'id'              => 340,
        'miseq_id'        => 21,
        'name'            => 'HUEDQ0501_B_KMT2C_1',
        'experiment_id'   => undef,
        'parent_plate_id' => undef,
        'gene'            => undef,
        'nhej_count'      => undef,
        'read_count'      => undef
    };
    my @exp_3_lost_wells = find_lost_wells( $test_exp_3, 'Miseq011' );

    my @expected_exp_3_lost_wells = ();
    is_deeply( \@exp_3_lost_wells, \@expected_exp_3_lost_wells,
        'No lost wells returned for test experiment 3' );

    return;
}

sub test_get_well_names : Test(3) {
    my $exp_1_wells          = get_well_names(165);
    my $expected_exp_1_wells = ['K09'];
    is_deeply( $exp_1_wells, $expected_exp_1_wells,
        'Correct well name retrieved for test experiment 1' );

    my $exp_2_wells          = ( get_well_names(167) );
    my @expected_exp_2_wells = ( 'A02', 'C09' );
    cmp_deeply(
        $exp_2_wells,
        bag(@expected_exp_2_wells),
        'Correct well names retrieved for test experiment 2'
    );

    my $exp_3_wells          = ( get_well_names(340) );
    my @expected_exp_3_wells = ( 'C02', 'F04', 'B04', '' );
    cmp_deeply(
        $exp_3_wells,
        bag(@expected_exp_3_wells),
        'Correct well names retrieved for test experiment 3'
    );

    return;
}

sub test_check_for_lost_well : Test(4) {
    my $lost_well_1 =
      check_for_lost_well( 'Miseq010', 195, 'Plate_2', ['K09'] );
    is( $lost_well_1, 'K01', 'Returns well name if well missing' );

    my $lost_well_2 =
      check_for_lost_well( 'Miseq010', 9, 'HUPEPD0017_1', [ 'C09', 'A02' ] );
    is( $lost_well_2, undef, 'Returns undef if well found' );

    my $lost_well_3 =
      check_for_lost_well( 'Miseq011', 30, 'HUEDQ0501_B_KMT2C_1',
        [ 'C02', 'F04', 'B04', '' ] );
    is( $lost_well_3, undef, 'Returns undef if well found' );

    my $lost_well_4 =
      check_for_lost_well( 'Miseq011', 2, 'DDD_SET1B_2_CAC', [] );
    is( $lost_well_4, undef, 'Returns undef if well and file missing' );

    return;
}

sub test_construct_file_path : Test(1) {
    my $file_path =
      construct_file_path( 'Test_plate_name', 123, 'Test_exp_name' );
    my $expected_file_path =
"$Bin/../../test_fixtures/quant_files/Test_plate_name/Test_exp_name/123.txt";
    is( $file_path, $expected_file_path, 'File path constructed correctly' );

    return;
}

sub test_create_report_file : Test(2) {
    my %test_mismatches = ( 'Exp' => ['A01'] );
    create_report_file( 'Test_plate', \%test_mismatches );
    my $test_plate_file = 'Test_plate_lost_wells.txt';
    ok( -e $test_plate_file, 'Report file created' );
    open( my $test_plate_fh, '<', $test_plate_file )
      or die "Couldn't open $test_plate_file: $!";
    chomp( my @test_plate_file_data = <$test_plate_fh> );
    close $test_plate_fh or die "Couldn't close $test_plate_file: $!";
    my @expected_test_plate_file_data = ('Exp: A01');
    is_deeply(
        \@test_plate_file_data,
        \@expected_test_plate_file_data,
        'Correct data in test report file'
    );

    return;
}

sub test_write_report : Test(2) {
    my %mismatches_1         = ();
    my $test_file_1_contents = '';
    open( my $test_fh_1, '>', $test_file_1_contents );
    write_report( $test_fh_1, \%mismatches_1 );
    close $test_fh_1;
    is( $test_file_1_contents, '', 'Report with no input data is empty' );

    my %mismatches_2 =
      ( 'Exp_1' => [ 'A01', 'B02', 'C03' ], 'Exp_2' => [ 'D04', 'E05' ] );
    my $test_file_2_contents = '';
    open( my $test_fh_2, '+<', \$test_file_2_contents );
    write_report( $test_fh_2, \%mismatches_2 );
    close $test_fh_2;
    ok(
        ( $test_file_2_contents eq "Exp_1: A01 B02 C03\nExp_2: D04 E05\n" )
          || (
            $test_file_2_contents eq "Exp_2: D04 E05\nExp_1: A01 B02 C03\n" ),
        'Correct data in test report'
    );

    return;
}

1;
