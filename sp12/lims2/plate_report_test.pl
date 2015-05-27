#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use LIMS2::Report;
use LIMS2::ReportGenerator::Plate;
use Log::Log4perl ':easy';
use Benchmark;
use Smart::Comments;

### [<time>] start

my $model = LIMS2::Model->new( user => 'tasks' );

Log::Log4perl->easy_init( { level => $INFO } );

my $plate_name = $ARGV[0];
my $plate = $model->retrieve_plate( { name => $plate_name } );
my $well_count = $plate->wells->count;
### $well_count
my $report_class = LIMS2::ReportGenerator::Plate->report_class_for( $plate->type_id );
$report_class =~ s/^.*\:\://;

my $t0 = Benchmark->new;

my $report_id = LIMS2::Report::generate_report(
    model      => $model,
    report     => $report_class,
    params     => { plate_id => $plate->id },
    async      => 0
);

die 'No report created' unless $report_id;

### $report_id

my ( $report_name, $report_fh ) = LIMS2::Report::read_report_from_disk( $report_id );

my $t1 = Benchmark->new;

my $td = timediff($t1,$t0);
print "the code took:",timestr($td),"\n";

### [<time>] name : $report_name
