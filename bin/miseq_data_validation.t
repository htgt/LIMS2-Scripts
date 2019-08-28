#!/usr/bin/env perl

use strict;
use warnings;

use Test::More tests => 5;

my @subs = qw(
                    miseq_data_validator
                    get_all_miseq_plates
                    get_miseq_exps 
                    find_lost_wells 
                    get_well_ids 
                    get_well_name 
                    check_for_lost_well 
                    construct_file_path 
                    create_report_file 
                    write_report 
);

BEGIN {
    use_ok( 'MiseqDataValidation', @subs );
}
