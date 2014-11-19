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
use DDP colored => 1;
use Const::Fast;
use YAML::Any qw( LoadFile );
use Math::Round qw( round );

use Smart::Comments;

=head2

1/ All protein coding genes for Mouse / Human.
2/ Get constitutive exons..
3/ Find target exons:
    a/ First 50% of coding region
    b/ exon must be larger than 200 bp ( no exon fragments of less that 100bp can be created )
4/ Within target exons grab search sequence
    a/ Leave 100bp of coding sequence either side
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

Output fill feed into next script that finds crisprs for insertions

FAILED TARGETS:
Info on genes we do not find intron insertion points for

=cut

my $log_level = $WARN;
GetOptions(
    'help'                 => sub { pod2usage( -verbose    => 1 ) },
    'man'                  => sub { pod2usage( -verbose    => 2 ) },
    'debug'                => sub { $log_level = $DEBUG },
    'verbose'              => sub { $log_level = $INFO },
    'trace'                => sub { $log_level = $TRACE },
    'genes-file=s'         => \my $genes_file,
    'gene=s'               => \my $single_gene,
    'species=s'            => \my $species,
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );
LOGDIE( 'Specify file with gene names' ) unless $genes_file;
LOGDIE( 'Must specify species' ) unless $species;

const my $DEFAULT_ASSEMBLY => $species eq 'Human' ? 'GRCh38' :  $species eq 'Mouse' ? 'GRCm38' : undef;
const my $DEFAULT_BUILD => 76;

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

# TODO add reason if I can
const my @FAILED_COLUMN_HEADERS => qw(
gene_id
marker_symbol
ensembl_id
ensembl_id_b
fail_reason
);

my $ensembl_util = WebAppCommon::Util::EnsEMBL->new( species => $species );
my $db = $ensembl_util->db_adaptor;
my $db_details = $db->to_hash;
WARN("Ensembl DB: " . $db_details->{'-DBNAME'});

my ( $output, $output_csv, $failed_output, $failed_output_csv );
$output = IO::File->new( 'targets.csv' , 'w' );
$output_csv = Text::CSV->new( { eol => "\n" } );
$output_csv->print( $output, \@COLUMN_HEADERS );

$failed_output = IO::File->new( 'failed.csv' , 'w' );
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
    my $ensembl_id = get_ensembl_id( $data );
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

    my @exons = @{ get_all_critical_exons( $gene, $data ) };
    unless ( @exons ) {
        ERROR('Unable to find any valid critical exons');
        $data->{fail_reason} = 'No critical exons found';
        push @failed_targets, $data;
        return;
    }

    # targets
    my $transcript = $gene->canonical_transcript;
    my $insertion_sites = get_intron_cassette_insertion_sites( \@exons, $transcript );
    unless ( %{ $insertion_sites } ) {
        ERROR('Unable to find any valid critical exons');
        $data->{fail_reason} = 'No insertion sites found';
        push @failed_targets, $data;
        return;
    }

    print_targets( $insertion_sites, $gene, $data );

    return;
}

sub get_ensembl_id {
    my $data = shift;

    if ( $data->{ensembl_id} ) {
        if ( $data->{ensembl_id_b} ) {
           if ( $data->{ensembl_id} eq $data->{ensembl_id_b} ) {
               return $data->{ensembl_id}
           }
           else {
               ERROR( 'Mismatch in ensembl ids: ' . $data->{ensembl_id} . ' and ' . $data->{ensembl_id_b});
               push @failed_targets, $data;
               return;
           }
        }
        else {
            return $data->{ensembl_id}
        }
    }
    else {
        if ( $data->{ensembl_id_b} ) {
            return $data->{ensembl_id_b}
        }
        else {
            ERROR( 'No Ensembl ID found' );
            push @failed_targets, $data;
            return;
        }
    }

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

=head2 get_all_critical_exons

All exons for the gene that are:
- constitutive ( belong to all coding transcipts )
- greater than 200 bases

Return list of exons sorted on ascending length

=cut
sub get_all_critical_exons {
    my ( $gene, $data ) = @_;

    return get_predefined_exons( $data->{exon_ids}, $gene, $data ) if $data->{exon_ids};

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

    my @ordered_critical_exons = sort most_five_prime @valid_critical_exons;

    return \@ordered_critical_exons;
}

=head2 get_predefined_exons

If exons targets have been pre-defined in the input use these
instead of trying to calculate the target exons.

=cut
sub get_predefined_exons {
    my ( $exon_ids, $gene, $data ) = @_;
    my @exons;
    my @exon_ids = split /,/, $data->{exon_ids};

    my %gene_exons = map{ $_->stable_id => 1 } @{ $gene->get_all_Exons };
    for my $exon_id ( @exon_ids ) {
        unless ( exists $gene_exons{$exon_id} ) {
            ERROR("The exon $exon_id can not be found on the gene: " . $gene->stable_id);
            next;
        }
        my $exon = $ensembl_util->exon_adaptor->fetch_by_stable_id( $exon_id );
        LOGDIE("Can not find specified exon $exon_id") unless $exon;
        push @exons, $exon;
    }

    return \@exons;
}

=head2 find_valid_exons

Exons that have greater than 200 coding bases, and coding.
Also must be within first 50% of coding sequence.

Create hash of valid exons, keyed on stable id
Create hash of exons for each transcript, keyed on transcript stable id

=cut
sub find_valid_exons {
    my ( $transcript, $valid_exons, $transcript_exons ) = @_;
    my @valid_exons;
    my @exons = @{ $transcript->get_all_Exons };

    my $cdna_mid = $transcript->cdna_coding_start + round( length($transcript->translateable_seq)/2 );

    for my $exon ( @exons ) {
        # skip exons which are non-coding
        unless ( $exon->coding_region_start( $transcript ) ) {
            DEBUG( 'Exon ' . $exon->stable_id . " is non coding in transcript " , $transcript->stable_id );
            next;
        }

        my $length = $exon->cdna_coding_end( $transcript ) - $exon->cdna_coding_start( $transcript ) + 1;
        if ( $length < 201 ) {
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

=head2 most_five_prime

Rank critical exons by following criteria by closest to 5 prime

=cut
## no critic(RequireFinalReturn)
sub most_five_prime {
    if ( $a->strand == 1 ) {
        $a->start <=> $b->start;
    }
    else {
        $b->start <=> $a->start;
    }
}
## use critic

=head2 get_intron_cassette_insertion_sites

4/ Within target exons grab search sequence
    a/ Leave 100bp of coding sequence either side
    b/ Cut off after reaching over 50% of coding region
- search seq: 100bp of exon seq to either side
    - AAGG / AAGA / CAGG / CAGA

        #my ( $start, $end ) = $transcript->cdna2genomic($cdna_mid, $cdna_mid);
=cut
sub get_intron_cassette_insertion_sites {
    my ( $exons, $transcript ) = @_;

    #get the mid point in cdna co-ordinates
    my $cdna_mid = $transcript->cdna_coding_start + round( length($transcript->translateable_seq)/2 );
    my ( $middle ) = $transcript->cdna2genomic($cdna_mid, $cdna_mid);
    $middle = $middle->start;

    my %valid_sites;
    my $strand = $transcript->seq_region_strand;

    for my $exon ( @{ $exons } ) {
        INFO( 'Searching for insertion sites in : ' . $exon->stable_id );

        if ( $exon->cdna_coding_start( $transcript ) > $cdna_mid ) {
            DEBUG( 'Exon ' . $exon->stable_id . " is starts more than 50% into coding bases" );
            next;
        }

        # must leave 100bp of coding sequence to either side of insertion site
        my $start = $exon->coding_region_start( $transcript ) + 100;
        my $end   = $exon->coding_region_end( $transcript ) - 100;
        # the two checks below should never occur as we have already filtered out
        # non-coding exons and small exons
        die ( 'Non coding exon' ) if !$start || !$end;
        die ( "Start $start greater than end $end" ) if $start > $end;

        if( $strand == 1 && $end > $middle ) {
            $end = $middle;
        }
        elsif ( $strand == -1 && $start < $middle  ) {
            $start = $middle;
        }

        # check we have enough sequence left to search for a insertion site
        if  ( ( $end - $start ) < 4 ) {
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

        my %sites;
        while ( $seq =~ /([AC]AG[GA])/g ) {
            my $ins_site = $1;
            my $end_pos = pos $seq;
            my $start_pos;
            if ( $strand == 1 ) {
                $start_pos = $start + ( $end_pos - 4 );
            }
            else {
                $start_pos = $start + ( $length - $end_pos );
            }
            $sites{$start_pos} = $ins_site;


        }

        if ( %sites ) {
            $valid_sites{$exon->stable_id}{length} = $exon->length;
            my ( $search_start, $search_end ) = minmax keys %sites;
            $search_start -= 1000;
            $search_end += 1000;
            my $slice = $ensembl_util->slice_adaptor->fetch_by_region(
                "chromosome",
                $transcript->seq_region_name,
                $search_start,
                $search_end,
            );
            my @sapI_sites = @{ map_SapI_binding_sites( $slice->seq, $search_start ) };

            for my $ins_site_start ( keys %sites ) {
                my $search_start = $ins_site_start - 1000;
                my $search_end = $ins_site_start + 1004;
                if ( any { $_ > $search_start && $_ < $search_end } @sapI_sites ) {
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
            my $num_sites = keys %{ $valid_sites{ $exon->stable_id }{sites} };
            INFO( ".. found $num_sites valid insertion sites" )
        }
        else {
            my $num_sites = keys %sites;
            INFO( ".. none found, have $num_sites insertion sites they all have flanking SapI sites" );
        }

    }

    return \%valid_sites;
}

=head2 map_SapI_binding_sites

Check for SapI binding sites within sequence
GCTCTTC
GAAGAGC

=cut
sub map_SapI_binding_sites {
    my ( $seq, $genomic_start ) = @_;

    my @binding_sites;
    while ( $seq =~ /(GCTCTTC|GAAGAGC)/g ) {
        my $bind_site = $1;
        my $end_pos = pos $seq;
        my $mid_pos = $genomic_start + ( $end_pos - 4 );
        #$binding_sites{$mid_pos} = $bind_site;
        push @binding_sites, $mid_pos;
    }
    #return \%binding_sites;
    return \@binding_sites;

    #my $c = () = $seq =~ /GCTCTTC|GAAGAGC/g;
    #return $c;

    #if ( index($seq, 'GCTCTTC') != -1 ) {
        #return 1;
    #}

    #if ( index($seq, 'GAAGAGC') != -1 ) {
        #return 1;
    #}

    return;
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

bk_design_targets.pl - Create design targets list given list of gene names.

=head1 SYNOPSIS

  gibson_design_targets.pl [options]

      --help            Display a brief help message
      --man             Display the manual page
      --debug           Debug output
      --verbose         Verbose output
      --trace           Trace output
      --genes-file      File with genes names.
      --gene            Specify only one gene from the file
      --species         Species of targets ( default Human )

      The genes file should be a csv file with 3 column headers: gene_id, marker_symbol and ensembl_id.
      The gene_id column will use HGNC ids or MGI ID's
      A fourth optionally column is exon_ids if the critical exons have been pre-defined.

=head1 DESCRIPTION


=head1 AUTHOR

Sajith Perera

=cut
