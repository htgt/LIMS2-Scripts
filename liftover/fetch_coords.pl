#!/usr/bin/perl
use strict;
use warnings;

use LIMS2::Model;
use IO::File;
use feature qw( say );
use Getopt::Long;

my $model = LIMS2::Model->new( user => 'tasks' );

my ( $species, $assembly );
GetOptions(
    'species=s'   => \$species,
    'assembly=s'  => \$assembly,
);

die ( 'Must specify --species' ) unless $species;
die ( 'Must specify --assembly' ) unless $assembly;

my @DATA = (
    {
        resultset => 'DesignOligoLocus',
        join      => { 'design_oligo' => 'design' },
        base      => 'design',
        id        => 'design_oligo_id',
    },
    {
        resultset => 'GenotypingPrimersLoci',
        join      => { 'genotyping_primer' => 'design' }, 
        base      => 'design',
        id        => 'genotyping_primer_id',
    },
    {
        resultset => 'CrisprLocus',
        join      => 'crispr', 
        base      => 'crispr',
        id        => 'crispr_id',
    },
);

for my $datum ( @DATA ) {
    say 'Working on: ' . $datum->{resultset};

    my @loci = $model->schema->resultset( $datum->{resultset} )->search(
        {
            $datum->{base} . '.species_id' => $species,
            assembly_id         => $assembly,
        },
        {
            join     => $datum->{join}, 
            prefetch => 'chr',
        }
    );

    print_bed_file(
        \@loci,
        $datum->{id},
        $datum->{resultset},
    );
}

say 'Working on: crispr primers'; 
my @all_crispr_primer_loci = $model->schema->resultset('CrisprPrimersLoci')->search(
    {},
    {
        join => { 'crispr_oligo' => [ 'crispr', 'crispr_pair', 'crispr_group' ] },
    }
);

my @crispr_primers;
for my $cpl ( @all_crispr_primer_loci ) {
    my $crispr;
    my $cp = $cpl->crispr_oligo;
    if ( $cp->crispr_id ) {
        $crispr = $cp->crispr;
    }
    elsif ( $cp->crispr_pair_id ) {
        $crispr = $cp->crispr_pair->left_crispr;
    }
    elsif ( $cp->crispr_group_id ) {
        $crispr = $cp->crispr_group->crisprs->first;
    }
    else {
        die( 'Crispr primer not linked to crisprs' );
    }

    push @crispr_primers, $cpl if $crispr->species_id eq $species;

}
print_bed_file( \@crispr_primers, 'crispr_oligo_id', 'CrisprPrimersLoci' ); 

sub print_bed_file {
    my ( $loci, $id, $name ) = @_;

    my $fh = IO::File->new( $name . '.bed' , 'w' );

    foreach my $locus ( @{ $loci } ) {
        say $fh join "\t",
            "chr" . $locus->chr->name,
            $locus->chr_start,
            $locus->chr_end,
            $locus->$id . ":" . $locus->chr_strand;
    }
}
