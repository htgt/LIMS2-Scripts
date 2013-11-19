#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Log::Log4perl ':easy';
use Pod::Usage;
use LIMS2::Model::Util::DesignInfo;
use Perl6::Slurp;
use feature qw( say );

my $log_level = $WARN;
GetOptions(
    'help'          => sub { pod2usage( -verbose => 1 ) },
    'man'           => sub { pod2usage( -verbose => 2 ) },
    'verbose'       => sub { $log_level = $INFO },
    'designs=s'     => \my $design_file,
    'plate=s'       => \my $plate_name,
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );
my $model = LIMS2::Model->new( user => 'lims2' );

my @designs;
if ( $plate_name ) {
    INFO( "Retrieving plate $plate_name" );
    my $plate = $model->retrieve_plate( { name => $plate_name } );

    @designs = map { $_->design } $plate->wells->all;
}
elsif ( $design_file ) {
    my @design_ids = map{ chomp; $_ } slurp $design_file;
    @designs = map{ $model->retrieve_design( { id => $_ } ) } @design_ids;
}
else {
    pod2usage(
        {   -message => 'Must specify either --designs or --plate',
            -verbose => 1,
        }
    );
}

say "gene,design_id,type,chromosome,strand,start,end,sequence";
for my $design ( @designs ) {
    design_report( $design );
}

sub design_report {
    my ( $design ) = @_;

    my $di = LIMS2::Model::Util::DesignInfo->new( { design => $design } );

    unless ( $di->type eq 'deletion' ) {
        ERROR( 'Unexpected design type: ' . $di->type );
        return;
    }
    my $genes = join( ' ', map{ $_->gene_id } $design->genes->all );
    my $oligos = $di->oligos;

    # +ve gene
    # U5 and G3 revcomp ( start / end swap ) 
    # -ve
    # D3 and G5 revcomp ( start / end swap )

    for my $oligo_type ( qw( G5 U5 D3 G3 ) ) {
        my $datum = $oligos->{$oligo_type};
        my ( $seq, $start, $end );
        if ( $di->chr_strand == 1 && $oligo_type =~ /U5|G3/ ) {
            $seq = revcomp( $datum->{seq} );
            $start = $datum->{end};
            $end = $datum->{start};
        }
        elsif ( $di->chr_strand == -1 && $oligo_type =~ /D3|G5/ ) {
            $seq = revcomp( $datum->{seq} );
            $start = $datum->{end};
            $end = $datum->{start};
        }
        else {
            $seq = $datum->{seq};
            $start = $datum->{start};
            $end = $datum->{end};
        }

        say $genes . ','
        . $design->id . ','
        . $oligo_type . ','
        . $datum->{chromosome} . ','
        . $datum->{strand} . ','
        . $start . ','
        . $end . ','
        . $seq;
    }
}

sub revcomp {
    my $dna     = shift;
    my $revcomp = reverse($dna);

    $revcomp =~ tr/ACGTacgt/TGCAtgca/;

    return $revcomp;
}


__END__

=head1 NAME

mirko_design_report.pl - basic design oligos report for Ross Cook

=head1 SYNOPSIS

  mirko_design_report.pl [options] plate-names

      --help            Display a brief help message
      --man             Display the manual page
      --verbose         Verbose output
      --designs         Specify file with list of designs
      --plate           Specify plate name

Send in a design plate or list of designs as arguments and the script will
produce a report listing the oligos of each design.

=head1 DESCRIPTION

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
