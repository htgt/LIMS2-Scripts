#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Log::Log4perl ':easy';
use LIMS2::Model::Util::CrisprESQC;

my $model = LIMS2::Model->new( user => 'lims2' );
Log::Log4perl->easy_init( { level => $INFO, layout => '%p %m%n' } );

my %params = (
    model                   => $model,
    plate_name              => 'HUEPD0005_1',
    sequencing_project_name => 'HUEPD0005_1_G',
    sub_seq_project         => 'HUEPD0005_1_G_1',
    commit                  => 0,
    user                    => 'sp12@sanger.ac.uk',
    species                 => 'Human',
);

my $qc_runner = LIMS2::Model::Util::CrisprESQC->new( %params );

my $qc_run = $qc_runner->analyse_plate;
