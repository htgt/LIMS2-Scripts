#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Test::MockModule;

use Test::More tests => 5;

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

my $module = Test::MockModule->new('MiseqDataValidation');
sub mock_construct_file_path {
    my ($plate_name, $i, $miseq_exp_name) = @_;
    my $file_path = "test_fixtures/quant_files/$plate_name/$miseq_exp_name/${i}.txt";
    return $file_path;
}
$module->mock('construct_file_path', &mock_construct_file_path); 

note('Testing package');

#BEGIN {
    use_ok( 'MiseqDataValidation', @subs );
#}

can_ok( __PACKAGE__, @subs );

note('Testing get_all_miseq_plates');

sub test_get_all_miseq_plates : Test(1) {
    my @test_plates = sort(get_all_miseq_plates());
    my @expected_plates = sort({'id' => 12307, 'name' => 'Miseq010', 'description' => '', 'type' => 'MISEQ', 'created_by' => '', 'created_at' => '', 'wells' => ['A02', 'C09', 'K09', 'O02']}, {'id' => 12316, 'name' => 'Miseq011', 'description' => '', 'type' => 'MISEQ', 'created_by' => '', 'created_at' => '', 'wells' => ['A01', 'B04', 'C02', 'C14', 'C17', 'E09', 'F04', 'I01', 'P09']});
    is_deeply( \@test_plates, \@expected_plates, 'Expected plate data retrieved');

    return;
}

note('Testing check_plate');

sub test_check_plate : Test(3) {
    my $test_miseq_10 = {'id' => 12307, 'name' => 'Miseq010', 'description' => '', 'type' => 'MISEQ', 'created_by' => '', 'created_at' => '', 'wells' => ['A02', 'C09', 'K09', 'O02']};
    check_plate($test_miseq_10);
    ok( -e 'Miseq010_lost_wells.txt', 'Report file created for Miseq 10' );

    my $line_count = 1;
    my %report_file_data;
    open(my $report_fh, '<', 'Miseq_10_lost_wells.txt') or die "Couldn't open Miseq 10 report file: $!";
    while (my $line = <$report_fh>) {
        $report_file_data{$line_count} = $line;
        $line_count++;
    }
    close $report_fh or die "Couldn't close Miseq 10 report file: $!";
    my %expected_file_data = (1 => 'CG_PDX_1_FPDM2: H05 M11', 2 => 'Plate_2: K01 O11', 3 => 'HUPEPD0017_1: M06');
    is_deeply(\%report_file_data, \%expected_file_data, 'Correct data in Miseq 10 report file');

    check_plate($test_miseq_11);
    ok( ! -e 'Miseq011_lost_wells.txt', 'Report file not created for Miseq 11' );

    return;
}

note('Testing get_mismatches');

sub test_get_mismatches : Test(2) {
    my $miseq_10_mismatches = get_mismatches(12307, 'Miseq010');
    my %expected_miseq_10_mismatches = ('CG_PDX_1_FPDM2' => ('H05', 'M11'), 'Plate_2' => ('K01', 'O11'), 'HUPEPD0017_1' => ('M06'));
    is_deeply( $miseq_10_mismatches, \%expected_miseq_10_mismatches, 'Correct mismatches for Miseq 10' );

    my $miseq_11_mismatches = get_mismatches(12316, 'Miseq011');
    my $expected_miseq_11_mismatches = {};
    is_deeply( $miseq_11_mismatches, $expected_miseq_11_mismatches, 'Correct mismatches for Miseq 11' );

    return;
}

note('Testing get_miseq_exps');

sub test_get_miseq_exps : Test(4) {
    my @miseq_10_exps = sort(get_miseq_exps(12307));
    my @expected_miseq_10_exps = sort({'id' => 170, 'miseq_id' => 13, 'name' => 'CG_PDX_1_FPDM2', 'experiment_id' => '', 'parent_plate_id' => '', 'gene' => '', 'nhej_count' => '', 'read_count' => ''}, {'id' => 165, 'miseq_id' => 13, 'name' => 'Plate_2', 'experiment_id' => '', 'parent_plate_id' => '', 'gene' => '', 'nhej_count' => '', 'read_count' => ''}, {'id' => 167, 'miseq_id' => 13, 'name' => 'HUPEPD0017_1', 'experiment_id' => '', 'parent_plate_id' => '', 'gene' => '', 'nhej_count' => '', 'read_count' => ''});  
    is_deeply( \@miseq_10_exps, \@expected_miseq_10_exps, 'Correct experiments for Miseq 10' ); 

    my @miseq_11_exps = sort(get_miseq_exps(12316));
    my @expected_miseq_11_exps = sort({'id' => 339, 'miseq_id' => 21, 'name' => 'DDD_SETD1B_2_CAC', 'experiment_id' => '', 'parent_plate_id' => '', 'gene' => '', 'nhej_count' => '', 'read_count' => ''}, {'id' => 340, 'miseq_id' => 21, 'name' => 'HUEDQ0501_B_KMT2C_1', 'experiment_id' => '', 'parent_plate_id' => '', 'gene' => '', 'nhej_count' => '', 'read_count' => ''}, {'id' => 288, 'miseq_id' => 21, 'name' => 'SETD5_2_PLATE11C', 'experiment_id' => '', 'parent_plate_id' => '', 'gene' => '', 'nhej_count' => '', 'read_count' => ''});  
    is_deeply( \@miseq_11_exps, \@expected_miseq_11_exps, 'Correct experiments for Miseq 11' ); 

    my @no_miseq_exps = get_miseq_exps(12345);
    my @expected_no_miseq_exps = ();
    is_deeply( \@no_miseq_exps, \@expected_no_miseq_exps, 'Returns empty list if no Miseq experiments linked to plate' );

    my @no_miseq_plate = get_miseq_exps(54321);
    my @expected_no_miseq_plate = ();
    is_deeply( \@no_miseq_plate, \@expected_no_miseq_plate, 'Returns empty list if no Miseq plate found' );

    return;
}

note ('Testing find_lost_wells');

sub test_find_lost_wells : Test(3) {
    my $test_exp_1 = {'id' => 165, 'miseq_id' => 13, 'name' => 'Plate_2', 'experiment_id' => '', 'parent_plate_id' => '', 'gene' => '', 'nhej_count' => '', 'read_count' => ''};
    my @exp_1_lost_wells = sort(find_lost_wells($test_exp_1, 'Miseq010'));
    my @expected_exp_1_lost_wells = ('K01', 'O11');
    is_deeply(\@exp_1_lost_wells, \@expected_exp_1_lost_wells, 'Correct lost wells returned for test experiment 1');

    my $test_exp_2 = {'id' => 167, 'miseq_id' => 13, 'name' => 'HUPEPD0017_1', 'experiment_id' => '', 'parent_plate_id' => '', 'gene' => '', 'nhej_count' => '', 'read_count' => ''};  
    my @exp_2_lost_wells = find_lost_wells($test_exp_2, 'Miseq010');
    my @expected_exp_2_lost_wells = ('M06');
    is_deeply(\@exp_2_lost_wells, \@expected_exp_2_lost_wells, 'Correct lost well returned for test experiment 2');

    my $test_exp_3 = {'id' => 340, 'miseq_id' => 21, 'name' => 'HUEDQ0501_B_KMT2C_1', 'experiment_id' => '', 'parent_plate_id' => '', 'gene' => '', 'nhej_count' => '', 'read_count' => ''};
    my @exp_3_lost_wells = find_lost_wells($test_exp_3, 'Miseq011');
    my @expected_exp_3_lost_wells = ();
    is_deeply(\@exp_3_lost_wells, \@expected_exp_3_lost_wells, 'No lost wells returned for test experiment 3');

    return;
}

note ('Testing get_well_names');

sub test_get_well_names : Test(3) {
    my @exp_1_wells = get_well_names(165);
    my @expected_exp_1_wells = ('K09');
    is_deeply(\@exp_1_wells, \@expected_exp_1_wells, 'Correct well name retrieved for test experiment 1');

    my @exp_2_wells = sort(get_well_names(167));
    my @expected_exp_2_wells = ('A02', 'C09');
    is_deeply(\@exp_2_wells, \@expected_exp_2_wells, 'Correct well names retrieved for test experiment 2');

    my @exp_3_wells = sort(get_well_names(340));
    my @expected_exp_3_wells = sort('C02', 'F04', 'B04', '');
    is_deeply(\@exp_3_wells, \@expected_exp_3_wells, 'Correct well names retrieved for test experiment 3');

    return;
}

note ('Testing check_for_lost_well');

sub test_check_for_lost_well : Test(3) {
    my $lost_well_1 = check_for_lost_well('Miseq010', 195, 'Plate_2', ['K09']);
    is($lost_well_1, 'K01', 'Returns well name if well missing');

    my $lost_well_2 = check_for_lost_well('Miseq010', 9, 'HUPEPD0017_1' ,['C09', 'A02']);
    is($lost_well_2, undef, 'Returns undef if well found');

    my $lost_well_3 = check_for_lost_well('Miseq011', 30, 'HUEDQ0501_B_KMT2C_1', ['C02', 'F04', 'B04', '']);
    is($lost_well_3, undef, 'Returns undef if well found');

    return;
}

note ('Testing construct_file_path');

sub test_contruct_file_path : Test(1) {
    my $file_path = contruct_file_path('Test_plate_name', 123, 'Test_exp_name');
    my $expected_file_path = 'test_fixtures/quant_files/Test_plate_name/Test_exp_name/123.txt'; 
    is($file_path, $expected_file_path, 'File path constructed correctly');
    return;
}

note ('Testing create_report_file');

sub test_create_report_file : Test(1) {
    my %test_mismatches = ('Exp_1' => ('A01', 'B02', 'C03'));
    create_report_file('Test_plate', \%test_mismatches);
    ok( -e 'Test_plate_lost_wells.txt', 'Report file created' );
    return;
}

note ('Testing write_report');

sub test_write_report : Test(2) {
    my $test_file_1 = 'Test_file.txt';
    my %mismatches = ();
    open(my $test_fh_1, '>', $test_file_1) or die "Couldn't open $test_file_1: $!";
    write_report($test_fh_1, \%mismatches);
    close $test_fh_1 or die "Couldn't close $test_file_1: $!";
    ok( -z $test_file_1, 'File with no input data is empty');

    my $test_file_2 = 'Test_plate_lost_wells.txt';
    my %test_file_2_data;
    my $line_count = 1;
    open(my $test_fh_2, '<', $test_file_2) or die "Couldn't open $test_file_2: $!";
    while (my $line = <$test_fh_2>) {
        $test_file_2_data{$line_count} = $line;
        $line_count++;
    }
    close $test_fh_2 or die "Couldn't close $test_file_2: $!";
    my %expected_test_file_2_data = (1 => 'Exp_1: A01 B02 C03');
    is_deeply(\%test_file_2_data, \%expected_test_file_2_data, 'Correct data in test report file');
    return;
}

unlink 'Test_file.txt', 'Test_plate_lost_wells.txt', 'Miseq010_lost_wells.txt', 'Miseq011_lost_wells.txt' or warn "Could not remove all files created by test: $!";
