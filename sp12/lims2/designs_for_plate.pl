#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Log::Log4perl ':easy';
use Pod::Usage;
use feature qw( say );

use Smart::Comments;

my $log_level = $WARN;
GetOptions(
    'help'          => sub { pod2usage( -verbose => 1 ) },
    'man'           => sub { pod2usage( -verbose => 2 ) },
    'debug'         => sub { $log_level = $DEBUG },
    'verbose'       => sub { $log_level = $INFO },
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );
my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

my $plate = $model->retrieve_plate( { name => $ARGV[0] } );

for my $well ( $plate->wells->all ) {
    my $design = $well->design;
    my $phase =  defined( $design->phase ) ? $design->phase : 'NONE';
    say "$well : " . $design->id . " - $phase";
}

__END__

=head1 NAME

designs_for_plate.pl - 

=head1 SYNOPSIS

  designs_for_plate.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --debug           Debug output
      --verbose         Verbose output

=head1 DESCRIPTION

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
