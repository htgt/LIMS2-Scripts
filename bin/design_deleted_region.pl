#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Pod::Usage;
use Perl6::Slurp;
use IO::Handle;
use WebAppCommon::Util::EnsEMBL;
use Bio::SeqIO;
use Bio::Seq;

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

GetOptions(
    'help'      => sub { pod2usage( -verbose => 1 ) },
    'man'       => sub { pod2usage( -verbose => 2 ) },
    'design=i'  => \my $design_id,
    'species=s' => \my $species_id,
);

die ( 'Must specify a species' ) unless $species_id;

my $ensembl_util = WebAppCommon::Util::EnsEMBL->new( species => $species_id );
my $sa = $ensembl_util->slice_adaptor;

my $output_fh = IO::Handle->new_from_fd( \*STDOUT, 'w' );
my $seq_out = Bio::SeqIO->new( -fh => $output_fh, -format => 'fasta' );

my @designs;
if ( $design_id ) {
    push @designs, $design_id;
}
elsif ( $ARGV[0] ) {
    push @designs, map{ chomp; $_ } slurp( $ARGV[0] );
}
else {
    podusage(2);
}

for my $id ( @designs ) {
    my $design = $model->c_retrieve_design( { id => $id } );

    die "unable to locate design $id" unless $design;
    die "Design $id is not $species_id" unless $design->species_id eq $species_id;

    my $start    = $design->info->cassette_start; 
    my $end      = $design->info->cassette_end; 
    my $chr_name = $design->info->chr_name; 

    my $slice = $sa->fetch_by_region(
        'chromosome',
        $chr_name,
        $start,
        $end,
    );

    my $seq = Bio::Seq->new(
        -display_id => $design->id,
        -alphabet   => 'dna',
        -seq        => $slice->seq,
    );
    $seq_out->write_seq( $seq );
}

=head1 NAME

delete_design.pl - create fasta file of designs deletion region 

=head1 SYNOPSIS

  delete_design.pl [options] [file]

      --help            Display a brief help message
      --man             Display the manual page
      --design          The id of the design
      --species         Species of the design ( must specify this )

  You must specify either a single design id through the command line
  or specify a file with a list of design ids, one per line.

=head1 DESCRIPTION

Grab the sequence for the deleted region for a design out output in a fasta file.
In this case the deleted region means:
    - Conditional: Between the U5 and U3 oligos ( non inclusive ).
    - Deletion: Between the U5 and D3 oligos ( non inclusive ).

=cut
