#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Try::Tiny;
use Text::CSV;
use Getopt::Long;
use feature qw(say);
use Bio::EnsEMBL::Registry;

#Setup Bio::EnsEMBL
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

my $slice_adaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'slice' );

my ( $file, $persist );
GetOptions(
    'file=s' => \$file,
    'persist'=> \$persist,
);

die( 'Specify file with design info' ) unless $file;
my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

open ( my $fh, '<', $file ) or die( "Can not open $file " . $! );
my $csv = Text::CSV->new();
$csv->column_names( @{ $csv->getline( $fh ) } );

while ( my $data = $csv->getline_hr( $fh ) ) {
    process_oligo( $data );
}

close $fh;

sub process_oligo {
    my $data = shift;

    my $design = $model->schema->resultset('Design')->find( { id => $data->{design_id} } );
    die('Can not find design ' . $data->{design_id}) unless $design;

    my $oligo = $design->oligos->find( { design_oligo_type_id => $data->{oligo_type} } );
    die('Can not find oligo ' . $data->{oligo_type} . ' ,design ' . $data->{design_id})
        unless $oligo;

    my $oligo_loci = $oligo->loci->find( { assembly_id => 'GRCm38' } );
    die('Can not find oligo loci on GRCm38 for ' . $data->{oligo_type} . ' ,design ' . $data->{design_id})
        unless $oligo_loci;

    my $chromosome = $model->schema->resultset('Chromosome')->find(
        {
            species_id => 'Mouse',
            name => $data->{chromosome},
        }
    );
    die('Can not find chromosome: ' . $data->{chromosome} )
        unless $chromosome;

    die( 'Start greater than end' ) if $data->{start} > $data->{end};
    die( 'Gap between start and end not 50' ) unless ( ( $data->{end} - $data->{start} ) + 1 ) == 50;

    my $sequence;
    if ( exists $data->{sequence} ) {
        die( 'Seq not 50 bases' ) unless length( $data->{sequence} ) == 50;
        $sequence = $data->{sequence};
    }
    else {
        my $slice = $slice_adaptor->fetch_by_region(
            'chromosome', $data->{chromosome}, $data->{start}, $data->{end} );
        $sequence = $slice->seq;
    }

    $model->txn_do(
        sub {
            try {
                $oligo->update( { seq => $sequence } );
                $oligo_loci->update(
                    {
                        chr_start  => $data->{start},
                        chr_end    => $data->{end},
                        chr_strand => $data->{strand},
                        chr_id     => $chromosome->id,
                    }
                );
                say 'Updated design oligo ' . $data->{oligo_type} . ' for design ' . $data->{design_id};
            }
            catch {
                say "Failed to update design $_ \n";
                $model->txn_rollback;
            };
            unless ( $persist ) {
                $model->txn_rollback;
                say 'Rollback!';
            }
        }
    );

    return;
}
