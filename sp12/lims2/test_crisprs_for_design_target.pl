#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use LIMS2::Model::Util::DesignTargets qw( bulk_crisprs_for_design_targets );
use YAML::Any;
use Smart::Comments;

my $model = LIMS2::Model->new( user => 'lims2' );


my $schema = $model->schema;

my $dt = $schema->resultset('DesignTarget')->find( { 'ensembl_exon_id' => $ARGV[0]  } );

return unless $dt;

my ( $crisprs, $crispr_pairs ) = bulk_crisprs_for_design_targets( $schema, [ $dt ], 'bwa' );

my $cps = $crispr_pairs->{ $dt->id };

my @cps = values %{ $cps };

my @good;
my @bad;
for my $cp ( @cps ) {
    if ( $cp->off_target_summary ) {
        my $summary = Load($cp->off_target_summary);
        if ( my $distance = $summary->{distance} ) {
            ### $distance
            $distance =~ s/\+//;
            if ( $distance >= 105 ) {
                push @good, $cp;
            }
            else {
                push @bad, $cp;
            }
        }
    }
    else {
        ### id no summary : $cp->id
    }

}
