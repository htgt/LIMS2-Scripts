#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Pod::Usage;
use Try::Tiny;
use JSON;
use Log::Log4perl ':easy';
use feature qw( say );

my $model = LIMS2::Model->new( user => 'lims2' );

GetOptions(
    'help'        => sub { pod2usage( -verbose => 1 ) },
    'man'         => sub { pod2usage( -verbose => 2 ) },
    'qc-run-id=s' => \my $qc_run_id,
    'all'         => \my $all,
);
Log::Log4perl->easy_init( { level => $INFO, layout => '%p %m%n' } );

my $qc_runs;
if ( $all ) {
    $qc_runs = $model->schema->resultset('CrisprEsQcRuns')->search(
        { species_id => 'Human' },
        { prefetch => {'crispr_es_qc_wells' => 'well'} }
    );
}
elsif ( $qc_run_id ) {
    $qc_runs = $model->schema->resultset('CrisprEsQcRuns')->search(
        { id => $qc_run_id },
        { prefetch => {'crispr_es_qc_wells' => 'well'} }
    );
}
else {
    die( 'Must specify qc-run-id or all' );
}


# PROBLEMS
# Not all variants called
# When no alignment then we can't call ( example: lots damage over crispr sites and reads stop )
# e.g. 8148C106-C0B3-11E4-8618-E1D03284BCA8 E04 / H07
# Some insertion / deletions in between crisprs not showing up ( not close enough to either cut site )
#

# CSV headers
say join ',', qw(
crispr_qc_well_id
well_name
gene
crispr_pair_id
chromosome
left_crispr_seq
left_crispr_start
left_crispr_end
left_crispr_cut_site
left_crispr_strand
right_crispr_seq
right_crispr_start
right_crispr_end
right_crispr_cut_site
right_crispr_strand
damage_type
var_start
var_type
var_ref_seq
var_alt_seq
var_length
var_read_depth
var_overlap_left_crispr
var_overlap_right_crispr
);


while (  my $qc_run = $qc_runs->next ) {
    INFO( 'Analysing run: ' . $qc_run->id );
    analyse_qc_run( $qc_run );
}

sub analyse_qc_run {
    my $qc_run = shift;

    my @report;
    for my $qc_well ( $qc_run->crispr_es_qc_wells ) {
        my $crispr = $qc_well->crispr;
        next unless $crispr;
        # only works on crispr pairs
        next unless $crispr->is_pair;

        my $right_crispr = $crispr->right_crispr;
        my $left_crispr = $crispr->left_crispr;
        my $right_crispr_locus = $right_crispr->current_locus;
        my $left_crispr_locus = $left_crispr->current_locus;

        # get formatted data for crispr es qc well
        my $well_data = $qc_well->format_well_data(
            sub { $model->find_genes( @_ ); }, #gene finder method
            {},
            $qc_run
        );

        # define cut sites for left and right crisprs.
        my $left_crispr_cut_site = $left_crispr_locus->chr_start + 5;
        my $right_crispr_cut_site = $right_crispr_locus->chr_end - 6;

        my @data = (
            $qc_well->id,                      # crispr qc well id
            $well_data->{well_name},           # well name
            $well_data->{gene},                # gene
            $crispr->id,                       # crispr pair id
            $left_crispr_locus->chr->name,     # chromosome
            $left_crispr->seq,                 # left crispr seq
            $left_crispr_locus->chr_start,     # left crispr start, PAM LEFT
            $left_crispr_locus->chr_end,       # left crispr end
            $left_crispr_cut_site,             # left crispr cut site
            $left_crispr_locus->chr_strand,    # left crispr strand
            $right_crispr->seq,                # right crispr seq
            $right_crispr_locus->chr_start,    # right crispr start, PAM RIGHT
            $right_crispr_locus->chr_end,      # right crispr end
            $right_crispr_cut_site,            # right crispr cut site
            $right_crispr_locus->chr_strand,   # right crispr strand
            $well_data->{damage_type},         # damage type
        );

        my $json = decode_json( $qc_well->analysis_data );
        my @vcf;

        # I only want the INDEL variants in the VCF file
        if ( exists $json->{non_merged_vcf} ) {
            @vcf = grep { /INDEL/ } grep{ !/#/ } split(/\n/, $json->{non_merged_vcf});
        }
        elsif ( $well_data->{has_vcf_file} ) {
            @vcf = grep { /INDEL/ } grep{ !/#/ } split(/\n/, $qc_well->vcf_file);
        }

        if ( @vcf ) {
            my $parsed_vcf = parse_vcf( \@vcf );
            for my $var ( @{ $parsed_vcf } ) {
                my @current_data = @data;
                push @current_data, $var->{start}, $var->{type}, $var->{ref}, $var->{alt}, $var->{length}, $var->{depth};
                push @current_data, var_overlap_cut_site( $var, $left_crispr_cut_site );
                push @current_data, var_overlap_cut_site( $var, $right_crispr_cut_site );
                push @report, \@current_data;
            }
        }
        else {
            push @data, 'No indel variants';
            push @report, \@data;
        }
    }

    # sort report on well name
    @report = sort { $a->[1] cmp $b->[1] } @report;
    say join( ',', map{ defined $_ ? $_ : '' } @{ $_ } ) for @report;
}


=head2 parse_vcf

Parse the INDEL variant lines from VCF file.
Work out length of INDEL.
Tab seperated data. 

=cut
sub parse_vcf {
    my $vcf = shift;

    my @variants;
    for my $var ( @{ $vcf } ) {
        my ( $chr, $start, $id, $ref, $alt, $qual, $filter, $info ) = split(/\t/, $var);
        my $type = length( $alt ) > length( $ref ) ? 'insertion' : 'deletion';
        my $var_length = abs( length( $ref ) - length( $alt ) );
        $info =~ /IDV=(?<depth>\d+)/;
        push @variants, {
            start  => $start,
            end    => $start + $var_length,
            ref    => $ref,
            alt    => $alt,
            type   => $type,
            length => $var_length,
            depth  => $+{depth},
        };
    }

    return \@variants;
}

=head2 var_overlaps_cut_site

Work out if variation overlaps the cut site.
It overlaps if it is within 3 bp of cut site.

=cut
sub var_overlap_cut_site {
    my ( $var, $cut_site ) = @_;
    if ( $var->{start} - 3 <= $cut_site && $var->{end} + 3 >= $cut_site ) {
        return 'yes';
    }
    return 'no';
}

__END__

=head1 NAME

crispr_vcf_damage_report.pl - Analyse damage over cut sites for crispr pairs

=head1 SYNOPSIS

  crispr_vcf_damage_report.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --qc-run-id       UUID of crispr es qc run
      --all             Run analysis against all crispr es qc runs

=head1 DESCRIPTION

Attempt to analyse damage caused by crispr pairs using VCF output from crispr es qc runs.
Figure out if there is damage overlapping the cut site for the left or right crispr.

Only works for crispr pairs at the moment.

=cut
