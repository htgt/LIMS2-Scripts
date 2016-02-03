#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Try::Tiny;
use Text::CSV;
use Getopt::Long;
use WebAppCommon::Util::EnsEMBL;
use Log::Log4perl ':easy';
use List::MoreUtils qw( all minmax any );
use IO::File;
use Pod::Usage;
use DDP colored => 0;
use Const::Fast;
use YAML::Any qw( LoadFile );
use Math::Round qw( round );

=head2

1/ All protein coding genes for Mouse / Human.
2/ Get constitutive exons ( also option to only get exons on canonical transcript )
3/ Find target exons:
    a/ First 50% of coding region
    b/ exon must be larger than 200 bp ( no exon fragments of less that 100bp can be created )
4/ Within target exons grab search sequence
    a/ Leave 80 / 100bp of coding sequence either side
    b/ Cut off after reaching over 50% of coding region
5/ Within this sequence search for insertion sites for intron cassette
    - AAGG / AAGA / CAGG / CAGA
    - A 1000bp either side of the insertion site should not have restriction enzyme site

TARGET OUTPUT:
BK wants to know number of genes that can be targeted with the above criteria
Input will also be feed into design creation, gibson designs, with inner primers flanking insertion point
- gene ( ensembl id, marker symbol, mgi id )
- exon ( exon id, rank, size )
- intron insertion ( sequence & coordinate )

Output will feed into next script ( bk_target_crispr_finder.pl ) that finds crisprs for insertions

FAILED TARGETS:
Info on genes we do not find intron insertion points for

=cut

my $log_level = $WARN;
GetOptions(
    'help'                 => sub { pod2usage( -verbose => 1 ) },
    'man'                  => sub { pod2usage( -verbose => 2 ) },
    'debug'                => sub { $log_level = $DEBUG },
    'verbose'              => sub { $log_level = $INFO },
    'trace'                => sub { $log_level = $TRACE },
    'genes-file=s'         => \my $genes_file,
    'gene=s'               => \my $single_gene,
    'species=s'            => \my $species,
    'canonical'            => \my $canonical,
    'coding-region=i'      => \my $custom_coding_region,
    'coding-flank=i'       => \my $custom_coding_flank,
    'sap1-flank=i'         => \my $custom_sap1_flank,
) or pod2usage(2);

my $CODING_REGION = $custom_coding_region ? ($custom_coding_region / 100) : 0.5;
my $CODING_FLANK = $custom_coding_flank || 60;
my $SAP1_FLANK = $custom_sap1_flank || 300;

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );
LOGDIE( 'Specify file with gene names' ) unless $genes_file;
LOGDIE( 'Must specify species' ) unless $species;

const my $DEFAULT_ASSEMBLY => $species eq 'Human' ? 'GRCh38' :  $species eq 'Mouse' ? 'GRCm38' : undef;
const my $DEFAULT_BUILD => 79;

LOGDIE( "Can not work out default assembly for species $species" ) unless $DEFAULT_ASSEMBLY;

WARN( "ASSEMBLY: $DEFAULT_ASSEMBLY, BUILD: $DEFAULT_BUILD" );

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
'ins_start',
'ins_end',
'ins_seq',
'assembly',
'build',
'species',
);

const my @FAILED_COLUMN_HEADERS => qw(
gene_id
marker_symbol
ensembl_id
fail_reason
critical_exons
);

my $ensembl_util = WebAppCommon::Util::EnsEMBL->new( species => $species );
my $db = $ensembl_util->db_adaptor;
my $db_details = $db->to_hash;
WARN("Ensembl DB: " . $db_details->{'-DBNAME'});

my ( $output, $output_csv, $failed_output, $failed_output_csv );
$output = IO::File->new( "targets_${genes_file}" , 'w' );
$output_csv = Text::CSV->new( { eol => "\n" } );
$output_csv->print( $output, \@COLUMN_HEADERS );

$failed_output = IO::File->new( "failed_targets_${genes_file}" , 'w' );
$failed_output_csv = Text::CSV->new( { eol => "\n" } );
$failed_output_csv->print( $failed_output, \@FAILED_COLUMN_HEADERS );

my @failed_targets;

## no critic(InputOutput::RequireBriefOpen)
{
    my $input_csv = Text::CSV->new();
    open ( my $input_fh, '<', $genes_file ) or die( "Can not open $genes_file " . $! );
    $input_csv->column_names( @{ $input_csv->getline( $input_fh ) } );
    while ( my $data = $input_csv->getline_hr( $input_fh ) ) {
        Log::Log4perl::NDC->remove;
        Log::Log4perl::NDC->push( $data->{gene_id} );
        Log::Log4perl::NDC->push( $data->{marker_symbol} );
        next if $single_gene && $single_gene ne $data->{marker_symbol};

        try{
            process_target( $data );
        }
        catch{
            ERROR('Problem processing target: ' . $_ );
            $data->{fail_reason} = $_;
            push @failed_targets, $data;
        };
    }
    for my $failed_target ( @failed_targets ) {
        $failed_output_csv->print( $failed_output, [ @{ $failed_target }{ @FAILED_COLUMN_HEADERS } ] );
    }
    close $input_fh;
}
## use critic

sub process_target {
    my $data = shift;
    my $ensembl_id = $data->{ensembl_id};
    return unless $ensembl_id;

    INFO( 'Target gene: ' . $ensembl_id );
    my $gene = $ensembl_util->gene_adaptor->fetch_by_stable_id( $ensembl_id );
    unless ( $gene ) {
        ERROR( "Can not find ensembl gene: " . $ensembl_id );
        $data->{fail_reason} = "Can not find gene $ensembl_id";
        push @failed_targets, $data;
        return;
    }

    # skip genes on patch / other non standard locations
    if ( $gene->seq_region_name !~ /^(?:\d+|X|Y|MT)$/ ) {
        ERROR( "Skipping gene on " . $gene->seq_region_name );
        $data->{fail_reason} = "Gene on chromosome: " . $gene->seq_region_name;
        push @failed_targets, $data;
        return;
    }

    my @exons;
    if ( $canonical ) {
        @exons = @{ get_valid_canonical_exons( $gene, $data ) };
    }
    else {
        @exons = @{ get_all_critical_exons( $gene, $data ) };
    }
    unless ( @exons ) {
        ERROR('Unable to find any valid exons');
        $data->{fail_reason} = 'No valid exons found';
        push @failed_targets, $data;
        return;
    }

    # targets
    my $transcript = $gene->canonical_transcript;
    my ( $insertion_sites, $site_count ) = get_intron_cassette_insertion_sites( \@exons, $transcript );
    unless ( %{ $insertion_sites } ) {
        ERROR('Unable to find any insertion sites');
        $data->{fail_reason} = "No valid insertion sites found ( $site_count )";
        $data->{critical_exons} = join( '|', map{ $_->stable_id } @exons );
        push @failed_targets, $data;
        return;
    }

    print_targets( $insertion_sites, $gene, $data );

    return;
}

=head2 print_targets

TARGET OUTPUT:
BK wants to know number of genes that can be targeted with the above criteria
Input will also be feed into design creation, gibson designs, with inner primers flanking insertion point
- gene ( ensembl id, marker symbol, mgi id )
- exon ( exon id, rank, size )
- intron insertion ( sequence & coordinate )

Output fill feed into next script that finds crisprs for insertions

=cut
sub print_targets {
    my ( $insertion_sites, $gene, $data ) = @_;

    my $canonical_transcript = $gene->canonical_transcript;
    my %base_params = (
        species              => $species,
        assembly             => $DEFAULT_ASSEMBLY,
        build                => $DEFAULT_BUILD,
        marker_symbol        => $data->{marker_symbol},
        gene_id              => $data->{gene_id},
        ensembl_gene_id      => $gene->stable_id,
        canonical_transcript => $canonical_transcript->stable_id,
    );

    for my $exon_id ( keys %{ $insertion_sites } ) {
        my %exon_params = %base_params;
        my $exon_rank = get_exon_rank( $exon_id, $canonical_transcript );

        $exon_params{ensembl_exon_id} = $exon_id;
        $exon_params{exon_size} = $insertion_sites->{$exon_id}{length};
        $exon_params{exon_rank} = $exon_rank;
        $exon_params{chr_name} = $canonical_transcript->seq_region_name;
        $exon_params{chr_strand} = $canonical_transcript->seq_region_strand;

        for my $ins_site ( keys %{ $insertion_sites->{$exon_id}{sites} } ) {
            my %target_params = %exon_params;
            $target_params{ins_start} = $ins_site;
            $target_params{ins_end}   = $ins_site + 3;
            $target_params{ins_seq}   = $insertion_sites->{$exon_id}{sites}{$ins_site};
            $output_csv->print( $output, [ @target_params{ @COLUMN_HEADERS } ] );
        }
    }

    return;
}

=head2 get_exon_rank

Get rank of exon on canonical transcript

=cut
sub get_exon_rank {
    my ( $exon_id, $canonical_transcript ) = @_;

    my $rank = 1;
    for my $current_exon ( @{ $canonical_transcript->get_all_Exons } ) {
        return $rank if $current_exon->stable_id eq $exon_id;
        $rank++;
    }

    return 0;
}

=head2 get_valid_canonical_exons


=cut
sub get_valid_canonical_exons {
    my ( $gene, $data ) = @_;

    my %valid_exons;
    my %transcript_exons;

    my $transcript = $gene->canonical_transcript;

    unless ( $transcript ) {
        WARN( 'Can not find canonical transcripts for gene: ' . $gene->stable_id );
        return [];
    }

    find_valid_exons( $transcript, \%valid_exons, \%transcript_exons );
    unless ( keys %valid_exons ) {
        WARN( 'No valid exons for gene: ' . $gene->stable_id );
        return [];
    }

    my @valid_exons = values %valid_exons;

    my $num_critical_exons = @valid_exons;
    INFO( "Has $num_critical_exons critical exons" );

    return \@valid_exons;
}

=head2 get_all_critical_exons

All exons for the gene that are:
- constitutive ( belong to all coding transcipts )
- greater than 200 bases

Return list of exons sorted on ascending length

=cut
sub get_all_critical_exons {
    my ( $gene, $data ) = @_;

    my %valid_exons;
    my %transcript_exons;
    my @coding_transcript_names;

    my @coding_transcripts = grep{ valid_coding_transcript($_) } @{ $gene->get_all_Transcripts };

    unless ( @coding_transcripts ) {
        WARN( 'Can not find coding transcripts for gene: ' . $gene->stable_id );
        return [];
    }

    for my $tran ( @coding_transcripts ) {
        push @coding_transcript_names, $tran->stable_id;
        find_valid_exons( $tran, \%valid_exons, \%transcript_exons );
    }
    unless ( keys %valid_exons ) {
        WARN( 'No valid exons for gene: ' . $gene->stable_id );
        return [];
    }
    DEBUG( 'Valid Exon Transcripts: ' . p( %transcript_exons ) );
    DEBUG( 'Valid Coding Transcripts: ' . p( @coding_transcript_names ) );

    my $critical_exons_ids = find_critical_exons( \%transcript_exons, \@coding_transcript_names );

    unless ( $critical_exons_ids ) {
        WARN( 'No critical exons for gene: ' . $gene->stable_id );
        return [];
    }

    my @valid_critical_exons;
    while ( my ( $exon_id, $exon ) = each %valid_exons ) {
        push @valid_critical_exons, $valid_exons{ $exon_id }
            if exists $critical_exons_ids->{ $exon_id };
    }

    unless ( @valid_critical_exons ){
        WARN( 'No valid exons that are also critical for gene ' . $gene->stable_id );
        return [];
    }

    my $num_critical_exons = @valid_critical_exons;
    INFO( "Has $num_critical_exons critical exons" );

    return \@valid_critical_exons;
}

=head2 find_valid_exons

Exons that have greater than 160 coding bases, and coding.
Also must be within first 50% of coding sequence.

Create hash of valid exons, keyed on stable id
Create hash of exons for each transcript, keyed on transcript stable id

=cut
sub find_valid_exons {
    my ( $transcript, $valid_exons, $transcript_exons ) = @_;
    my @valid_exons;
    my @exons = @{ $transcript->get_all_Exons };

    # if cannot find start of coding, exit
    return unless $transcript->cdna_coding_start;

    my $cdna_mid = $transcript->cdna_coding_start + round( $CODING_REGION * length($transcript->translateable_seq) );

    for my $exon ( @exons ) {
        # skip exons which are non-coding
        unless ( $exon->coding_region_start( $transcript ) ) {
            DEBUG( 'Exon ' . $exon->stable_id . " is non coding in transcript " , $transcript->stable_id );
            next;
        }

        my $length = $exon->cdna_coding_end( $transcript ) - $exon->cdna_coding_start( $transcript ) + 1;
        if ( $length < (2 * $CODING_FLANK) ) {
            DEBUG( 'Exon ' . $exon->stable_id . " only has $length coding bases " );
            next;
        }

        # filter out exons that are wholly in the last 50% of coding bases
        # a more specific check is done later for exons that straddle the 50% mark
        if ( $exon->cdna_coding_start( $transcript ) > $cdna_mid ) {
            DEBUG( 'Exon ' . $exon->stable_id . " is starts more than 50% into coding bases" );
            next;
        }

        TRACE( 'Exon ' . $exon->stable_id . ' is VALID' );
        push @valid_exons, $exon;
    }

    for my $exon ( @valid_exons ) {
        push @{ $transcript_exons->{ $exon->stable_id } }, $transcript->stable_id;
        $valid_exons->{ $exon->stable_id } = $exon
            unless exists $valid_exons->{ $exon->stable_id };
    }

    return;
}

=head2 find_critical_exons

Find exons that belong to all coding transcripts

=cut
sub find_critical_exons {
    my ( $transcript_exons, $coding_transcript_names ) = @_;
    my %critical_exons;

    for my $exon ( keys %{ $transcript_exons } ) {
        my %transcripts = map{ $_ => 1  } @{ $transcript_exons->{$exon} };

        if ( all { exists $transcripts{ $_ } } @{ $coding_transcript_names } ) {
            TRACE( "Exon $exon is critical" );
            $critical_exons{ $exon } = 1;
        }
        else {
            TRACE( "Exon $exon not present in all valid transcripts" );
        }
    }
    return unless keys %critical_exons;

    return \%critical_exons;
}

=head2 valid_coding_transcript

Return true if the transcript is a valid transcript by our measure:
Not nonsense mediated decay.
CDS region in complete ( has proper start and end ).

=cut
sub valid_coding_transcript {
    my ( $transcript ) = @_;
    my $id = $transcript->stable_id;

    TRACE( "$id biotype: " . $transcript->biotype );
    if ( !$transcript->translation ) {
        TRACE("Transcript $id is non protein coding");
        return 0;
    }

    if ( $transcript->biotype eq 'nonsense_mediated_decay') {
        TRACE("Transcript $id is NMD");
        return 0;
    }

    # CDS incomplete check, both 5' and 3'
    if ( _get_transcript_attribute( $transcript, 'cds_end_NF' ) ) {
        TRACE("Transcript $id has incomplete CDS end");
        return 0;
    }

    if ( _get_transcript_attribute( $transcript, 'cds_start_NF' ) ) {
        TRACE("Transcript $id has incomplete CDS start");
        return 0;
    }

    TRACE("Transcript $id is VALID");
    return 1;
}

=head2 get_intron_cassette_insertion_sites

We have a array of critical exons that we can search through now.

Search for the intron cassette insertion sites within each exon:
    a/ Leave 100bp of coding sequence either side of insertion site
    b/ Insertion site must be within first 50% of coding bases
    c/ Insertion sites: AAGG / AAGA / CAGG / CAGA

Once sites found make sure there are not any SapI binding sites with 1000 bases
to either side.

=cut
sub get_intron_cassette_insertion_sites {
    my ( $exons, $transcript ) = @_;

    my $count_sites = 0;
    #get the mid point in cdna co-ordinates
    my $cdna_mid = $transcript->cdna_coding_start + round( length($transcript->translateable_seq)/2 );
    my ( $middle ) = $transcript->cdna2genomic($cdna_mid, $cdna_mid);
    $middle = $middle->start;

    my %valid_sites;
    my $strand = $transcript->seq_region_strand;

    for my $exon ( @{ $exons } ) {
        INFO( 'Searching for insertion sites in : ' . $exon->stable_id );

        # must leave at least coding sequence to either side of insertion site
        my $start = $exon->coding_region_start( $transcript ) + $CODING_FLANK;
        my $end   = $exon->coding_region_end( $transcript ) - $CODING_FLANK;

        # insertion site must be within first 50% of coding sequences
        if( $strand == 1 && $end > $middle ) {
            $end = $middle;
        }
        elsif ( $strand == -1 && $start < $middle  ) {
            $start = $middle;
        }

        # check we have enough sequence left to search for a insertion site
        if ( ( $end - $start ) < 4 ) {
            DEBUG( 'Coding region in exon not big enough for insertion site' );
            next;
        }

        # need sequence on correct strand
        my $slice = $ensembl_util->slice_adaptor->fetch_by_region(
            "chromosome",
            $transcript->seq_region_name,
            $start,
            $end,
            $strand,
        );
        my $seq = $slice->seq;
        my $length = length( $seq );

        # search for intron cassette insertion sites in sequence
        my %sites;
        while ( $seq =~ /([AC]?AG[GA])/g ) {
            my $ins_site = $1;
            my $end_pos = pos $seq;

            my $start_pos;
            if ( $strand == 1 ) {
                $start_pos = $start + ( $end_pos - length($ins_site) );
            }
            else {
                $start_pos = $start + ( $length - $end_pos );
            }
            $sites{$start_pos} = $ins_site;
            $count_sites++;
        }

        # If we have insertion sites filter out ones that have SapI sites within 1000 bases
        # either side
        if ( %sites ) {
            $valid_sites{$exon->stable_id}{length} = $exon->length;
            my ( $search_start, $search_end ) = minmax keys %sites;
            $search_start -= $SAP1_FLANK;
            $search_end   += $SAP1_FLANK;
            my $site_slice = $ensembl_util->slice_adaptor->fetch_by_region(
                "chromosome",
                $transcript->seq_region_name,
                $search_start,
                $search_end,
            );
            # get array of SapI binding site coordinates
            my @sapI_sites = @{ map_SapI_binding_sites( $site_slice->seq, $search_start ) };

            for my $ins_site_start ( keys %sites ) {
                my $region_start = $ins_site_start - $SAP1_FLANK;
                my $region_end   = $ins_site_start + $SAP1_FLANK + length($sites{$ins_site_start}); # the cassette insertion site size

                if ( any { $_ > $region_start && $_ < $region_end } @sapI_sites ) {
                    next;
                }
                $valid_sites{$exon->stable_id}{sites}{$ins_site_start} = $sites{ $ins_site_start };
            }
        }
        else {
            INFO( '.. none found' );
            next;
        }

        if ( exists $valid_sites{ $exon->stable_id }{sites} ) {
            my $num_valid_sites = keys %{ $valid_sites{ $exon->stable_id }{sites} };
            INFO( ".. found $num_valid_sites valid insertion sites" )
        }
        else {
            my $num_sites = keys %sites;
            INFO( ".. none found, have $num_sites insertion sites they all have flanking SapI sites" );
            delete $valid_sites{ $exon->stable_id };
        }

    }

    return ( \%valid_sites, $count_sites);
}

=head2 map_SapI_binding_sites

Check for SapI binding sites within sequence
GCTCTTC
GAAGAGC

Return array of coordinates for the binding sites, if any

=cut
sub map_SapI_binding_sites {
    my ( $seq, $genomic_start ) = @_;

    my @binding_sites;
    while ( $seq =~ /(GCTCTTC|GAAGAGC)/g ) {
        my $bind_site = $1;
        my $end_pos = pos $seq;
        my $mid_pos = $genomic_start + ( $end_pos - 4 );
        push @binding_sites, $mid_pos;
    }

    return \@binding_sites;
}

sub _get_transcript_attribute {
    my ( $transcript, $code ) = @_;

    my ( $attr ) = @{ $transcript->get_all_Attributes($code) };
    if ( $attr ) {
        return $attr->value();
    }
    return 0;
}

__END__

=head1 NAME

bk_target_finder.pl - Find custom target sites for BK  ( Bon-Kyoung )

=head1 SYNOPSIS

  bk_target_finder.pl [options]
      --help            Display a brief help message
      --man             Display the manual page
      --debug           Debug output
      --verbose         Verbose output
      --trace           Trace output
      --genes-file      File with genes names.
      --gene            Specify only one gene from the file
      --species         Species of targets ( default Human )
      --canonical       Use exons from canonical transcript instead of critical exons
      --coding-region   Exon must be in first coding-region percentage (default 50)
      --coding-flank    Number of coding bases in exon that must flank insertion site (default 60)
      --sap1-flank      Number of coding bases flanking insertion site that must not contain a Sap1 site (default 300)

      The genes file should be a csv file with 3 column headers: gene_id, marker_symbol and ensembl_id.
      The gene_id column will use HGNC ids or MGI ID's
      A fourth optionally column is exon_ids if the critical exons have been pre-defined.

=head1 DESCRIPTION


=cut
