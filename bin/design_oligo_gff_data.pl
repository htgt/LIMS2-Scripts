#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Perl6::Slurp;
use Log::Log4perl ':easy';
use Getopt::Long;
use Const::Fast;
use feature qw( say );

my $log_level = $WARN;

GetOptions(
    debug   => sub { $log_level = $DEBUG },
    verbose => sub { $log_level = $INFO },
) and @ARGV == 1
    or die "Usage: $0 FILE\n";

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

const my @FIELDS => qw(
seqname
source
feature
start
end
score
strand
frame
group
);

my %designs = map { split /,/ } map { chomp; $_ } slurp( $ARGV[0] );

print_gff_line( { map{ $_ => ucfirst($_) } @FIELDS } );

for my $gene ( keys %designs ) {
    my $design_id = $designs{$gene};
    Log::Log4perl::NDC->push( $design_id );
    INFO( "Generate gff data for design $design_id, gene $gene" );

    my $design = $model->retrieve_design( { id => $design_id } );
    unless( $design ) {
        ERROR( "Unable to find design $design_id" );
        next;
    }

    generate_gff_output( $design, $gene );

    Log::Log4perl::NDC->remove;
}

sub generate_gff_output {
    my ( $design, $gene ) = @_;

    my $default_locus = $design->default_locus;

    my %design_data = (
        seqname => 'chr' . $default_locus->chr_name,
        source  => $gene,
        feature => 'target',
        group   => $gene,
    );

    print_gff_line( \%design_data );

    for my $oligo ( qw( g5 u5 d3 g3 ) ) {
        generate_oligo_gff_output( $oligo, $default_locus );
    }
}

sub generate_oligo_gff_output {
    my ( $oligo, $default_locus ) = @_;

    my $start_name = $oligo . '_start';
    my $end_name = $oligo . '_end';
    my %oligo_data = (
        seqname => 'chr' . $default_locus->chr_name,
        source  => uc( $oligo ),
        feature => 'oligo',
        start   => $default_locus->$start_name,
        end     => $default_locus->$end_name,
        group   => $default_locus->design_id,
    );

    print_gff_line( \%oligo_data );
}

sub print_gff_line {
    my $data = shift;

    my $line;
    for my $field ( @FIELDS ) {
        if ( exists $data->{$field} ) {
            $line .= $data->{$field};
        }
        else {
            $line .= '.';
        }
        $line .= "\t";
    }

    say $line;
}

=head1 SYNOPSIS

design_oligo_gff_data.pl[options] gene_design_list.csv

      --verbose         Print informational messages
      --debug           Print debug messages

=head1 DESCRIPTION

Expects a csv file, 2 values for each row, first is gene name, second design_id.

Outputs gff data to STDOUT.

NOTE: Only works for deletion / insertion designs currently

=head1 AUTHOR

Sajith Perera

=head1 BUGS

=headl TODO

Make it work for conditional designs

=cut

