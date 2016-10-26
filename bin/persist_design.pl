#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Getopt::Long;
use Pod::Usage;
use LIMS2::REST::Client;
use YAML::Any qw( LoadFile DumpFile );
use Moose;
use Data::Dumper;
use Log::Log4perl qw( :easy );

my ( $file, $persist, $dir );
GetOptions(
    'help'            => sub { pod2usage( -verbose => 1 ) },
    'man'             => sub { pod2usage( -verbose => 2 ) },
    'file=s'          => \$file,
    'debug'           => \my $debug,
) or pod2usage(2);

die( 'Specify file with design info' ) unless $file;

has lims2_api => (
    is         => 'ro',
    isa        => 'LIMS2::REST::Client',
    traits     => [ 'NoGetopt' ],
    lazy_build => 1
);

sub _build_lims2_api {
    my $self = shift;

    return LIMS2::REST::Client->new_with_config();
}

my $lims = {
    lims2_api         => _build_lims2_api(),
    dir               => $dir,
    design_method     => 'conditional-inversion',
};

my $design_data = YAML::Any::LoadFile( $file );
my $self = $lims->{lims2_api}->_build_ua();
my $design = $lims->{lims2_api}->POST('design', $design_data );
print ('Design persisted: ' . $design->{id} );

__END__

=head1 NAME

create-multiple-designs.pl - Create multiple designs

=head1 SYNOPSIS

  create-multiple-designs.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --debug           Print debug messages
      --file            File with design details.
      --persist         Persist newly created designs to LIMS2
      --alt-designs     Create alternate designs
      --dir             Directory where design-create output goes
      --gene            Only create this gene(s), picked from input file
      --param           Specify additional param(s) not in file
      --dry-run         Just print out command that would be called, don't call it
      --param           Specify additional parameter(s) to send to design creation program

=head1 DESCRIPTION

Takes design information for multiple designs from file and tries to create these designs.

The file will be a csv file, each row represents a design, each column a parameter given to the
design create command. The column headers represent the parameter name.

=cut
