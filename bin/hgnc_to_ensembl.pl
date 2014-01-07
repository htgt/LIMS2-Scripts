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
    'help'    => sub { pod2usage( -verbose => 1 ) },
    'man'     => sub { pod2usage( -verbose => 2 ) },
    'debug'   => sub { $log_level = $DEBUG },
    'verbose' => sub { $log_level = $INFO },
    'file=s'  => \my $file,
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %m%n' } );

LOGDIE( 'Specify file with gene names' ) unless $file;

my @hgnc_gene_names = map{ chomp; $_ } slurp($file);

my $ensembl_util = LIMS2::Util::EnsEMBL->new( species => 'Human' );

for my $hgnc_name ( @hgnc_gene_names ) {
    INFO("Gene: $hgnc_name");

    my $gene = get_ensembl_gene( $hgnc_name );
    my $output = "$hgnc_name,";
    $output .= $gene->stable_id if $gene;

    say $output;
}

sub get_ensembl_gene {
    my $hgnc_name = shift;
    my $ensembl_gene;

    my @genes = @{ $ensembl_util->gene_adaptor->fetch_all_by_external_name($hgnc_name, 'HGNC') };

    unless( @genes ) {
        WARN( "Unable to find gene $hgnc_name in EnsEMBL" );
        return;
    }

    if ( scalar(@genes) > 1 ) {
        DEBUG("Found multiple EnsEMBL genes for $hgnc_name");

        my @filtered_genes = grep{ $_->external_name eq $hgnc_name } @genes;
        for my $gene ( @genes ) {
            DEBUG($gene->external_name );
            DEBUG($gene->stable_id );

        }
        unless( @filtered_genes ) {
            WARN('None of the genes have same HGNC name');
            return;
        }

        if ( scalar(@filtered_genes) == 1 ) {
            DEBUG("Only one has correct HGNC symbol");
            $ensembl_gene = shift @filtered_genes;
        }
        else {
            DEBUG("Can not identify single EnsEMBL gene linked to $hgnc_name");
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

hgnc_to_ensembl.pl - Convert HGNC ids to EnsEMBL ids

=head1 SYNOPSIS

  hgnc_to_ensembl.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --debug           Debug output
      --verbose         Verbose output
      --file            File with genes names.

=head1 DESCRIPTION

Input newline seperated file of HGCN symbols, output csv of those
HGNC symbols linked with corresponding EnsEMBL ids;

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
