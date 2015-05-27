#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Benchmark;
use Smart::Comments;

### [<time>] start

my $model = LIMS2::Model->new( user => 'tasks' );

my $plate_name = $ARGV[0];
my $plate = $model->retrieve_plate( { name => $plate_name } );

my $t0 = Benchmark->new;

my $rs = $model->schema->resultset( 'PlateReport' )->search(
    {},
    {
        prefetch => 'well', 
        bind => [ $plate->id ],
    }
);

my $data = $rs->consolidate;

my $t1 = Benchmark->new;

my $td = timediff($t1,$t0);
print "the code took:",timestr($td),"\n";

