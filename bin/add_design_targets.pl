#!/usr/bin/env perl

use strict;
use warnings;

use feature qw( say );

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Try::Tiny;

my ( @exon_ids, $species, $assembly, $comment, $build );
GetOptions(
    'help'            => sub { pod2usage( 1 ) },
    'man'             => sub { pod2usage( 2 ) },
    'exon-ids=s{,}'   => \@exon_ids,
    'species=s'       => sub { my ( $name, $val ) = @_; $species = ucfirst(lc $val); },
    'assembly=s'      => \$assembly,
    'comment=s'       => \$comment,
    'build=s'         => \$build,

) or pod2usage( 2 );

pod2usage( 1 ) unless @exon_ids;

use LIMS2::Model;
use LIMS2::Util::EnsEMBL;

my $model = LIMS2::Model->new( user => 'lims2' );
my $ens = LIMS2::Util::EnsEMBL->new( species => $species );

#check species is in lims2 db
die "Couldn't find '$species' in species db" 
    unless $model->schema->resultset('Species')->find( { id => $species } );


#see if any design targets already exist
my @design_targets = $model->schema->resultset("DesignTarget")->search(
    { 
        ensembl_exon_id => { -IN => [ @exon_ids ] },
        build_id        => $build,
        assembly_id     => $assembly
    }
);

die "The following targets already exist: " . join(",", map { $_->ensembl_exon_id } @design_targets) 
    if @design_targets;

my %id_map;
for my $exon_id ( @exon_ids ) {
    say "Processing $exon_id";
    my $gene = $ens->gene_adaptor->fetch_by_exon_stable_id( $exon_id );

    die "Couldn't find a gene for $exon_id!" unless $gene;

    my $rank = get_rank( $gene->canonical_transcript, $exon_id );
    my $exon = $ens->exon_adaptor->fetch_by_stable_id( $exon_id );
    my $chr_id = $model->schema->resultset("Chromosome")->find( {
        species_id => $species,
        name    => $exon->seq_region_name
    } )->id;

    my $data = {
        automatically_picked => 1,
        canonical_transcript => $gene->canonical_transcript->stable_id,
        marker_symbol   => $gene->external_name,
        ensembl_gene_id => $gene->stable_id,
        species_id      => $species,
        assembly_id     => $assembly,
        build_id        => $build,
        ensembl_exon_id => $exon_id,
        exon_size       => $exon->seq->length,
        exon_rank       => $rank,
        chr_id          => $chr_id,
        chr_start       => $exon->seq_region_start,
        chr_end         => $exon->seq_region_end,
        chr_strand      => $exon->seq_region_strand,
        comment         => $comment,
    };

    my $dt = $model->schema->resultset("DesignTarget")->create( $data );

    $id_map{$exon_id} = $dt->id;
}

#display all the ids
say "Exon id - dt id";
say $_ . " - " . $id_map{$_} for keys %id_map;

#should change one in util ensembl to do this
sub get_rank {
    my ( $transcript, $exon_id ) = @_;

    my $rank = 1;
    for my $exon ( @{ $transcript->get_all_Exons } ) {
        if ( $exon->stable_id eq $exon_id ) {
            return $rank;
        }

        $rank++;
    }

    say "Couldn't find $exon_id in " . $transcript->stable_id;
    
    #return undef if not found
    return;
}

1;

__END__

=head1 NAME

add_design_targets.pl - add design targets for a list of exons

=head1 SYNOPSIS

add_design_targets.pl [options]
               
    --exon-ids         the exons you want to add design targets for
    --species          the species these exons belong to
    --assembly         the assembly the design target will belong to
    --build            the ensembl build
    --comment          additional information for the exon(s) [optional]
    --help             show this dialog

Example usage:

perl add_design_targets.pl --species human --build 73 --assembly GRCh37 --exon-ids ENSE00001353200 ENSE00001140991

=head1 DESCRIPTION

Given an exon id (or list of exon ids) this script will delete every crispr associated (+/- 200bp), and
the related tables.

=head AUTHOR

Alex Hodgkins

=cut