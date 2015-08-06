#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use WGE::Model::DB;
use Const::Fast;
use Data::Dumper;
use IO::File;
use Getopt::Long;
use Text::CSV;
use Log::Log4perl ':easy';
use YAML::Any;
use Pod::Usage;
use Try::Tiny;

my $log_level = $WARN;
GetOptions(
    'help'          => sub { pod2usage( -verbose => 1 ) },
    'man'           => sub { pod2usage( -verbose => 2 ) },
    'debug'         => sub { $log_level = $DEBUG },
    'verbose'       => sub { $log_level = $INFO },
    'trace'         => sub { $log_level = $TRACE },
    'target-file=s' => \my $targets_file,
    'gene=s'        => \my $single_gene,
    'species=s'     => \my $species,
    'flanking=i'    => \my $flanking,
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );

LOGDIE( 'Specify file with targets' ) unless $targets_file;
LOGDIE( 'Must specify species' ) unless $species;
$species = 'Grch38' if $species eq 'Human';

my $wge = WGE::Model::DB->new;
INFO( "Finding species $species in WGE");
my $SPECIES_ID = $wge->resultset('Species')->find( { id => $species } )->numerical_id;
die "Couldn't find species $species" unless $SPECIES_ID;

my $input_csv = Text::CSV->new();
open ( my $input_fh, '<', $targets_file ) or die( "Can not open $targets_file " . $! );
$input_csv->column_names( @{ $input_csv->getline( $input_fh ) } );

my @failed_targets;
my @original_fields = $input_csv->column_names;
push @original_fields, 'fail_reason';
my $failed_output = IO::File->new( 'failed_crispr_targets.csv' , 'w' );
my $failed_output_csv = Text::CSV->new( { eol => "\n" } );
$failed_output_csv->print( $failed_output, \@original_fields );

const my @COLUMN_HEADERS => (
'gene_id',
'marker_symbol',
'ensembl_gene_id',
'ensembl_exon_id',
'exon_size',
'exon_rank',
'canonical_transcript',
'chr_name',
'chr_strand',
'insertion_seq',
'insertion_start',
'insertion_end',
'assembly',
'build',
'species',
'id',
'chr_name',
'chr_start',
'chr_end',
'seq',
'pam_right',
'off_target_summary',
);

my $output = IO::File->new( 'crisprs_for_targets.csv' , 'w' );
my $output_csv = Text::CSV->new( { eol => "\n" } );
$output_csv->print( $output, \@COLUMN_HEADERS );

{
    while ( my $data = $input_csv->getline_hr( $input_fh ) ) {
        Log::Log4perl::NDC->remove;
        Log::Log4perl::NDC->push( $data->{marker_symbol} );
        Log::Log4perl::NDC->push( $data->{ensembl_exon_id} );
        next if $single_gene && $single_gene ne $data->{marker_symbol};

        try{
            find_target_crisprs( $data );
        }
        catch{
            ERROR('Problem processing target: ' . $_ );
            $data->{fail_reason} = $_;
            push @failed_targets, $data;
        };
    }

    for my $failed_target ( @failed_targets ) {
        $failed_output_csv->print( $failed_output, [ @{ $failed_target }{ @original_fields } ] );
    }
    close $input_fh;
}

=head2 find_target_crisprs

desc

=cut
sub find_target_crisprs {
    my ( $data ) = @_;

    my @crisprs;
    if ( $flanking ) {
        @crisprs = crisprs_flanking_region(
            {
                chr_name  => $data->{chr_name},
                chr_start => $data->{ins_start},
                chr_end   => $data->{ins_end},
                strand    => $data->{chr_strand},
                seq       => $data->{ins_seq},
            }
        );
    }
    else {
        @crisprs = crisprs_overlapping_region(
            {
                chr_name  => $data->{chr_name},
                chr_start => $data->{ins_start},
                chr_end   => $data->{ins_end},
                strand    => $data->{chr_strand},
                seq       => $data->{ins_seq},
            }
        );
    }

    unless ( @crisprs ) {
        $data->{fail_reason} = 'No crisprs found overlapping insertion site';
        push @failed_targets, $data;
        WARN( '.. no crisprs found' );
        return;
    }

    my @valid_crisprs = grep{ is_valid_crispr( $_ ) } map{ $_->as_hash } @crisprs;
    unless ( @valid_crisprs ) {
        my $num_crisprs = @crisprs;
        $data->{fail_reason} = "Found $num_crisprs crispr(s) but non are valid";
        push @failed_targets, $data;
        WARN( ".. found $num_crisprs crisprs but non valid" );
        return;
    }
    my $num_valid_crisprs = @valid_crisprs;
    DEBUG( ".. found $num_valid_crisprs valid crisprs" );

    print_target( $data, \@valid_crisprs );
}

=head2 crisprs_overlapping_region

Diagram This...
Numbers allow insertion site to cut crispr site from between 2nd / 3rd bases and the 14th / 15th ones.

=cut
sub crisprs_overlapping_region {
    my ( $params ) = @_;

    DEBUG("Getting crisprs overlapping: " . region_str( $params ) );

    my ( $pam_right_start_offset, $pam_right_end_offset, $pam_left_start_offset, $pam_left_end_offset );
    if ( $params->{strand} == 1 ) {
        # FOR +ve stranded genes
        $pam_right_start_offset = 18; # last 2 bases can be GG, which should not match PAM in crispr
        $pam_right_end_offset   = 6;
        $pam_left_start_offset  = 11;
        $pam_left_end_offset    = -1; # first 2 bases in ins seq not CC so can leave at 0
    }
    else {
        # FOR -ve stranded genes
        $pam_right_start_offset = 20; # last 2 bases can be GG, which should not match PAM in crispr
        $pam_right_end_offset   = 8;
        $pam_left_start_offset  = 13;
        $pam_left_end_offset    = 1; # first 2 bases in ins seq not CC so can leave at 0
    }

    my @pam_right_crisprs = $wge->resultset('Crispr')->search(
        {
            'species_id'  => $SPECIES_ID,
            'chr_name'    => $params->{chr_name},
            'pam_right'   => 1,
            'chr_start'   => {
                -between => [ $params->{chr_start} - $pam_right_start_offset, $params->{chr_start} - $pam_right_end_offset ],
            },
        },
    );

    my @pam_left_crisprs = $wge->resultset('Crispr')->search(
        {
            'species_id'  => $SPECIES_ID,
            'chr_name'    => $params->{chr_name},
            'pam_right'   => 0,
            'chr_start'   => {
                -between => [ $params->{chr_start} - $pam_left_start_offset, $params->{chr_start} - $pam_left_end_offset ],
            },
        },
    );

    return ( @pam_right_crisprs, @pam_left_crisprs );
}

=head2 crisprs_flanking_region


=cut
sub crisprs_flanking_region {
    my ( $params ) = @_;

    DEBUG("Getting crisprs flanking $flanking bases " . region_str( $params ) );

    my @crisprs = $wge->resultset('Crispr')->search(
        {
            'species_id'  => $SPECIES_ID,
            'chr_name'    => $params->{chr_name},
            'chr_start'   => {
                -between => [ $params->{chr_start} - $flanking, $params->{chr_end} + $flanking ],
            },
        },
    );

    return @crisprs;
}

sub region_str {
    my $crispr = shift;

    unless ( $crispr->{chr_name} && $crispr->{chr_start} && $crispr->{chr_end} ) {
        die Dumper( $crispr );
    }

    return $crispr->{chr_name} . ':' . $crispr->{chr_start} . '-' . $crispr->{chr_end};
}

=head2 is_valid_crispr

desc

=cut
sub is_valid_crispr {
    my ( $crispr ) = @_;

    # TODO filter crisprs
    die "off_target_summary is null for " . Dumper( $crispr ) unless $crispr->{off_target_summary};

    # TODO exact matches only - check this
    return Load( $crispr->{off_target_summary} )->{0} <= 1;
}

=head2 print_target


=cut
sub print_target {
    my ( $data, $crisprs ) = @_;

    my %base_params = (
        species              => $species,
        marker_symbol        => $data->{marker_symbol},
        gene_id              => $data->{gene_id},
        ensembl_gene_id      => $data->{ensembl_gene_id},
        canonical_transcript => $data->{canonical_transcript},
        ensembl_exon_id      => $data->{ensembl_exon_id},
        chr_name             => $data->{chr_name},
        chr_strand           => $data->{chr_strand},
        exon_size            => $data->{exon_size},
        exon_rank            => $data->{exon_rank},
        insertion_seq        => $data->{ins_seq},
        insertion_start      => $data->{ins_start},
        insertion_end        => $data->{ins_end},
        assembly             => $data->{assembly},
        build                => $data->{build},
    );

    for my $crispr ( @{ $crisprs } ) {
        my %crispr_params = %base_params;

        $crispr_params{id}                 = $crispr->{id};
        $crispr_params{chr_name}           = $crispr->{chr_name};
        $crispr_params{chr_start}          = $crispr->{chr_start};
        $crispr_params{chr_end}            = $crispr->{chr_end};
        $crispr_params{seq}                = $crispr->{seq};
        $crispr_params{pam_right}          = $crispr->{pam_right};
        $crispr_params{off_target_summary} = $crispr->{off_target_summary};

        $output_csv->print( $output, [ @crispr_params{ @COLUMN_HEADERS } ] );
    }

    return;
}
