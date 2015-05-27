#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Try::Tiny;
use YAML::Any;
use Getopt::Long;
use feature qw(say);
use Bio::EnsEMBL::Registry;

use Smart::Comments;

#Setup Bio::EnsEMBL
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

my $slice_adaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'slice' );

my ( $file, $persist );
GetOptions(
    'persist'=> \$persist,
);

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

for my $file ( @ARGV ) {
    my $data = YAML::Any::LoadFile( $file );

    say 'WORKING ON: ' . $data->{id};
    update_design_oligos( $data );
}


sub update_design_oligos {
    my $data = shift;

    my $design = $model->schema->resultset('Design')->find( { id => $data->{id} } );
    die('Can not find design ' . $data->{id}) unless $design;

    for my $oligo_data ( @{ $data->{oligos} } ) {
        update_oligo( $oligo_data, $design );
    }
}

sub update_oligo {
    my ( $oligo_data, $design ) = @_;

    my $oligo = $design->oligos->find( { design_oligo_type_id => $oligo_data->{type} } );
    die('Can not find oligo ' . $oligo_data->{type} . ' ,design ' . $design->id)
        unless $oligo;

    my $loci_data = shift @{ $oligo_data->{loci} };

    my $oligo_loci = $oligo->loci->find( { assembly_id => 'GRCm38' } );
    die('Can not find oligo loci on GRCm38 for ' . $oligo_data->{type} . ' ,design ' . $design->id)
        unless $oligo_loci;

    my $chromosome = $model->schema->resultset('Chromosome')->find(
        {
            species_id => 'Mouse',
            name => $loci_data->{chr_name},
        }
    );
    die('Can not find chromosome: ' . $loci_data->{chr_name} )
        unless $chromosome;

    die( 'Start greater than end' ) if $loci_data->{chr_start} > $loci_data->{chr_end};
    die( 'Gap between start and end not 50' )
        unless ( ( $loci_data->{chr_end} - $loci_data->{chr_start} ) + 1 ) == 50;

    die( 'Seq not 50 bases' )
        unless length( $oligo_data->{seq} ) == 50;

    $model->txn_do(
        sub {
            try {
                $oligo->update( { seq => $oligo_data->{seq} } );
                $oligo_loci->update(
                    {
                        chr_start  => $loci_data->{chr_start},
                        chr_end    => $loci_data->{chr_end},
                        chr_strand => $loci_data->{chr_strand},
                        chr_id     => $chromosome->id,
                    }
                );
                say 'Updated design oligo ' . $oligo_data->{type} . ' for design ' . $design->id;
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
