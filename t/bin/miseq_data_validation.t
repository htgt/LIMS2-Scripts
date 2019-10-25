#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';

BEGIN {
    use Log::Log4perl qw( :easy );
    Log::Log4perl->easy_init($FATAL);
}

use FindBin qw($Bin);
use lib "$Bin/../lib";
use Test::Class;
use MiseqDataValidationTest;

Test::Class->runtests;

unlink 'Test_plate_lost_wells.txt', 'Miseq010_lost_wells.txt',
  'Miseq011_lost_wells.txt'
  or warn "Could not remove all files created by test: $!";

1;
