#!/usr/bin/perl
use strict;
use warnings;

use LIMS2::Model;
use Getopt::Long;
use feature qw( say );
use Try::Tiny;

my %data = (
    DesignOligoLocus      => 'design_oligo_id',
    GenotypingPrimersLoci => 'genotyping_primer_id',
    CrisprLocus           => 'crispr_id',
    CrisprPrimersLoci     => 'crispr_oligo_id',
);

my ( $species, $assembly, $resultset, $file, $commit );
GetOptions(
    'species=s'   => \$species,
    'assembly=s'  => \$assembly,
    'resultset=s' => \$resultset,
    'data-file=s' => \$file,
    'commit'      => \$commit,
);

$assembly //= 'GRCh38';
$species  //= 'Human';

die ( 'Must specify --resultset' ) unless $resultset;
die ( "Unknown resultset: $resultset" ) unless exists $data{ $resultset };

my $model = LIMS2::Model->new( user => 'tasks' );

my %chr_ids;
my @chromosomes = $model->schema->resultset("Chromosome")->search( { species_id => $species } );
for my $chr ( @chromosomes ) {
    $chr_ids{ 'chr' . $chr->name } = $chr->id;
}

$model->txn_do(
    sub {
        try{
            open (my $fh, "<", $file) or die $!;
            update_loci( $fh );
            unless ( $commit ) {
                print "non-commit mode, rollback\n";
                $model->txn_rollback;
            }
        }
        catch {
            print "failed: $_\n";
            $model->txn_rollback;
        };
    }
);

sub update_loci {
    my ( $fh ) = shift;

    my $count = 0;
    foreach my $line ( <$fh> ){
        $count++;
        chomp $line;
        my ($chr,$start,$end, $identifier) = split "\t", $line;
        my ($oligo_id, $strand) = split ":", $identifier;

        my $params = {
            $data{ $resultset } => $oligo_id,
            assembly_id         => $assembly,
            chr_start           => $start,
            chr_end             => $end,
            chr_strand          => $strand,
            chr_id              => $chr_ids{ $chr },
        };

        say "Adding new $resultset: $oligo_id ( $count )";
        my $locus = $model->schema->resultset( $resultset )->create($params);

        unless ( $locus->chr_start == $start and $locus->chr_end == $end ) {
            die "Start and end coords not as expected. Expected: $start, $end, got: "
                . $locus->chr_start . ", "
                . $locus->chr_end . "\n";
        }
    }
}
