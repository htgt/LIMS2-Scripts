#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Perl6::Slurp;
use Try::Tiny;
use Getopt::Long;
use LIMS2::Util::EnsEMBL;
use Log::Log4perl ':easy';
use Pod::Usage;
use feature qw( say );

my $log_level = $WARN;
GetOptions(
    'help'      => sub { pod2usage( -verbose => 1 ) },
    'man'       => sub { pod2usage( -verbose => 2 ) },
    'debug'     => sub { $log_level = $DEBUG },
    'verbose'   => sub { $log_level = $INFO },
    'species=s' => \my $species,
    'type=s'    => \my $gene_type,
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );

my $ensembl_util = LIMS2::Util::EnsEMBL->new( species => $species );

while ( my $gene_name = <> ) {
    chomp($gene_name);
    Log::Log4perl::NDC->remove;
    Log::Log4perl::NDC->push( $gene_name );

    my $gene = get_ensembl_gene( $gene_name, $gene_type );
    my $output = "$gene_name,";
    $output .= $gene->stable_id . ',' if $gene;
    $output .= $gene->external_name if $gene;

    say $output;
}

sub get_ensembl_gene {
    my ( $gene_name, $type ) = @_;
    my $ensembl_gene;

    my $id = $gene_name;
    if ( $type eq 'HGNC' ) {
        $gene_name =~ /HGNC:(\d+)/;
        $id = $1;
    }

    my @genes = @{ $ensembl_util->gene_adaptor->fetch_all_by_external_name($id, $type) };

    unless( @genes ) {
        WARN( 'Unable to find gene EnsEMBL' );
        return;
    }

    if ( scalar(@genes) > 1 ) {
        DEBUG('Found multiple EnsEMBL genes');

        my @filtered_genes = grep{ $_->external_name eq $gene_name } @genes;
        for my $gene ( @genes ) {
            DEBUG($gene->stable_id );
            DEBUG($gene->external_name );

            #my @dbentries = @{ $gene->get_all_DBEntries( 'MGI' ) };
            #for my $d ( @dbentries ) {
                #say "$_ : " . $d->$_ for qw( primary_id optional_id dbname display_id description dbID description ensembl_object_type  );
            #}

        }
        unless( @filtered_genes ) {
            WARN("None of the genes have same $type name");
            return;
        }

        if ( scalar(@filtered_genes) == 1 ) {
            DEBUG("Only one has correct $type symbol");
            $ensembl_gene = shift @filtered_genes;
        }
        else {
            DEBUG("Can not identify single linked EnsEMBL gene");
            $ensembl_gene = pick_from_multiple_genes( \@filtered_genes );
        }
    }
    else {
        $ensembl_gene = shift @genes;
    }
    INFO('EnsEMBL gene: ' . $ensembl_gene->stable_id );

    return $ensembl_gene;
}

sub pick_from_multiple_genes {
    my $genes = shift;

    for my $g ( @{ $genes } ) {
        if ( $g->stable_id =~ /^ENSG\d+$/ ) {
            return $g;
        }
    }
    WARN('Can not pick the right EnsEMBL gene');
}

__END__

=head1 NAME

hgnc_to_ensembl.pl - Convert external gene ids to EnsEMBL ids

=head1 SYNOPSIS

  hgnc_to_ensembl.pl [options] <gene names>

      --help            Display a brief help message
      --man             Display the manual page
      --debug           Debug output
      --verbose         Verbose output
      --species         Species for gene, Human or Mouse
      --type            Type of external gene names you are specifying ( e.g MGI, HGNC )

=head1 DESCRIPTION

Input newline seperated list of gene names.
Outputs a csv of those gene names linked with corresponding EnsEMBL ids.

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
