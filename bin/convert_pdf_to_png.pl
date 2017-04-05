#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Getopt::Long;
use Pod::Usage;
use PostScript::Convert;
use Cwd;
use File::Find;

#my ( $filename );
#GetOptions(
#    'file=s'          => \$filename,
#) or pod2usage(2);
$DB::single=1;
my $base = getcwd;
my $dirs = [];
my $wanted = sub { _wanted($dirs) };

find($wanted, $base);
foreach my $file (@$dirs) {
    psconvert($file, format => 'png');
    print ("Converted: " . $file . "\n");
}
sub _wanted {
   return if ! -e; 
   my ($dirs) = @_;

   push( @$dirs, $File::Find::name ) if $File::Find::name=~ /1b\.Indel_size_distribution_percentage\.pdf/;
}
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
