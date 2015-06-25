#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Getopt::Long;
use Log::Log4perl ':easy';
use YAML::Any;
use Perl6::Slurp;
use List::Util qw( first );
use LIMS2::Model::Util::DesignInfo;
use LIMS2::Model::Util::DesignTargets qw(
design_target_report_for_genes
);

use LIMS2::Model;

use Smart::Comments;

my $model = LIMS2::Model->new( user => 'tasks' );

my $log_level = $WARN;
GetOptions(
    'debug'  => sub { $log_level = $DEBUG },
    'gene=s' => \my @input_genes,
);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );

my $schema = $model->schema;

#my @genes;

#if ( @input_genes ) {
    #@genes = @input_genes
#}
#elsif ( my $input_file = $ARGV[0] ) {
    #@genes = map{ chomp; $_ } slurp $ARGV[0];
#}
#else {
    #@genes = ( 'HGNC:5157' , 'HGNC:21222', 'HGNC:21427' );
#}

#my $gene_string = join("\n", @genes);
#my ( $report_data, $sorted_genes ) = design_target_report_for_genes( $model->schema, $gene_string, 'Human' );


#my @dts = $schema->resultset( 'DesignTarget' )->search( { marker_symbol => 'CDH1', species_id => 'Human' } );

#my @data = $schema->resultset( 'DesignTargetCrisprs' )->search( {
        #design_target_id => { 'IN' => [ map { $_->id } @dts ] },
    #}
#);

#my %exon_crispr_count;
#for my $d ( @data ) {
    #$exon_crispr_count{ $d->ensembl_exon_id }++;
#}

## %exon_crispr_count

my $design = $schema->resultset('Design')->find( {
        id => 1001884
    },
    {
        prefetch => { 'oligos' => { 'loci' => 'chr' } }
    }
);


my %design_oligos_data;
for my $oligo ( $design->oligos ) {
    my ( $locus ) = grep{ $_->assembly_id eq 'GRCh37' } $oligo->loci;

    my %oligo_data = (
        start      => $locus->chr_start,
        end        => $locus->chr_end,
        chromosome => $locus->chr->name,
        strand     => $locus->chr_strand,
    );
    $oligo_data{seq} = $oligo->seq;

    $design_oligos_data{ $oligo->design_oligo_type_id } = \%oligo_data;
}

my $di = LIMS2::Model::Util::DesignInfo->new(
    design => $design,
    default_assembly => 'GRCh37',
    oligos => \%design_oligos_data,
);

my $oligo_hash = $di->oligos;

### $oligo_hash

### ts : $di->target_region_start
### te : $di->target_region_end


__END__

=head1 NAME

design_crispr_data.pl -

=head1 SYNOPSIS

  design_crispr_data.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --debug           Debug output
      --verbose         Verbose output
      --trace           Trace output

=head1 DESCRIPTION

=head1 AUTHOR

Sajith Perera

=head1 BUGS

None reported... yet.

=head1 TODO

=cut
