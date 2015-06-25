#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use EngSeqBuilder;
use YAML::Any qw( LoadFile );
use Log::Log4perl ':easy';

use Smart::Comments;

my $esb = EngSeqBuilder->new(
    configfile            => $ENV{ENG_SEQ_BUILDER_CONF},
    max_vector_seq_length => 250000
);
Log::Log4perl->easy_init( { level => $TRACE, layout => '%p %x %m%n' } );

my $params = LoadFile( $ARGV[0] );

for my $well ( values %{ $params->{wells} } ) {

    my $method = $well->{eng_seq_method};
    my $params = $well->{eng_seq_params};

    delete $params->{target_region_start};
    delete $params->{target_region_end};

    ### $method
    ### $params

    my $gbk = $esb->$method( $params );
    last;
}
