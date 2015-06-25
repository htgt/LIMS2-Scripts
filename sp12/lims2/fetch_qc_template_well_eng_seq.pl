#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Log::Log4perl ':easy';
use Pod::Usage;
use IO::File;
use IO::Handle;
use LIMS2::Model::Util::EngSeqParams qw( generate_genbank_for_qc_well );

my $log_level = $WARN;
GetOptions(
    'help'          => sub { pod2usage( -verbose => 1 ) },
    'man'           => sub { pod2usage( -verbose => 2 ) },
    'debug'         => sub { $log_level = $DEBUG },
    'verbose'       => sub { $log_level = $INFO },
    'well=s'        => \my $well,
    'plate=s'       => \my $plate,
) or pod2usage(2);

die( "Must specify both --well and --plate" ) if !$well || !$plate;

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );
my $model = LIMS2::Model->new( user => 'lims2' );

my $qc_well = $model->schema->resultset( 'QcTemplateWell' )->find(
    {
        'qc_template.name' => $plate,
        name               => $well,
    },
    {
        join => 'qc_template',
    }
);

my $ofh = IO::Handle->new->fdopen( fileno(STDOUT), 'w' );
generate_genbank_for_qc_well( $qc_well, $ofh );

__END__

=head1 NAME

test_crispr_template.pl - generate qc template well genbank file 

=head1 SYNOPSIS

  test_crispr_template.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --debug           Debug output
      --verbose         Verbose output
      --well            Name of qc template well
      --plate           Name of qc template

=head1 DESCRIPTION

Generate genbank file for a given qc template well.
Must specify both a plate name and a well name.

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
