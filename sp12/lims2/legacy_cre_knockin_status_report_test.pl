#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use LIMS2::Model::Util::LegacyCreKnockInProjectReport;
use Benchmark;
use Smart::Comments;
use Log::Log4perl ':easy';

Log::Log4perl->easy_init( { level => $INFO } );

my $model = LIMS2::Model->new( user => 'tasks' );

my $t0 = Benchmark->new;

my $cre_ki_report_generator = LIMS2::Model::Util::LegacyCreKnockInProjectReport->new(
    model      => $model,
    project_id => 579,
);
$cre_ki_report_generator->generate_report_data;
my $report_data = $cre_ki_report_generator->report_data;

### $report_data

my $t1 = Benchmark->new;

my $td = timediff($t1,$t0);
#print "the code took:",timestr($td),"\n";

