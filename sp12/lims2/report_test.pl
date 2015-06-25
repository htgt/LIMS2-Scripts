#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use LIMS2::ReportGenerator::Plate;
use Log::Log4perl ':easy';

use Benchmark;
use Smart::Comments;

my $model = LIMS2::Model->new( user => 'tasks' );

Log::Log4perl->easy_init( { level => $INFO } );

my $plate = $model->schema->resultset('Plate')->find( { id => 2613 } );

my $p0 = Benchmark->new;

my $report_class = LIMS2::ReportGenerator::Plate->report_class_for( $plate->type_id );
### $report_class

my $report = $report_class->new( model => $model, species => 'Mouse', plate => $plate );

my $p1 = Benchmark->new;
my $pt = timestr( timediff($p1,$p0) );
### Time : $pt

