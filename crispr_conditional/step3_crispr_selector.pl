#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use WGE::Model::DB;
use Const::Fast;
# use Data::Dumper;
# use IO::File;
use Getopt::Long;
use Text::CSV;
# use Log::Log4perl ':easy';
# use YAML::Any;
# use Pod::Usage;
use Try::Tiny;


use Bio::EnsEMBL::Registry;

use Smart::Comments;



my $species = 'Grch38';

my $wge = WGE::Model::DB->new;
# INFO( "Finding species $species in WGE");
my $SPECIES_ID = $wge->resultset('Species')->find( { id => $species } )->numerical_id;
die "Couldn't find species $species" unless $SPECIES_ID;

$species = 'Human' if $species eq 'Grch38';



# csv file with a list of targets (inserion sites) with primers
my $targets_file = shift @ARGV;


my $input_csv = Text::CSV->new();
open ( my $input_fh, '<', $targets_file ) or die( "Can not open $targets_file " . $! );
$input_csv->column_names( @{ $input_csv->getline( $input_fh ) } );

# my @original_fields = $input_csv->column_names;
# ## @original_fields


const my @COLUMN_HEADERS => (
'gene_id',
'marker_symbol',
'assembly',
'build',
'species',
'ensembl_gene_id',
'ensembl_exon_id',
'exon_size',
'exon_rank',
'canonical_transcript',
'exon_start',
'exon_end',
'chr_name',
'chr_strand',
'ins_start',
'ins_end',
'ins_seq',
'insertion_site',
'split_exon_L',
'split_exon_R',
'LAL_seq',
'LAL_loc',
'LAL_gc',
'LAL_tm',
'LAR_seq',
'LAR_loc',
'LAR_gc',
'LAR_tm',
'LA_size',
'RAL_seq',
'RAL_loc',
'RAL_gc',
'RAL_tm',
'RAR_seq',
'RAR_loc',
'RAR_gc',
'RAR_tm',
'RA_size',

'crispr_id',
'crispr_start',
'crispr_end',
'crispr_seq',
'pam_right',
'proximity',
'off_target_summary',
'off_target_count',

);



# my $output = IO::File->new( "crisprs_${targets_file}" , 'w' );


open(my $output_fh, ">", "crisprs_${targets_file}");

my $output_csv = Text::CSV->new( { eol => "\n" } );

$output_csv->print( $output_fh, \@COLUMN_HEADERS );

close($output_fh);

# get ensembl ready to query
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous');

my $data;

while ( $data = $input_csv->getline_hr( $input_fh ) ) {

    # my $ensembl_util = WebAppCommon::Util::EnsEMBL->new( species => $species );
    # my $db = $ensembl_util->db_adaptor;
    # my $db_details = $db->to_hash;
    # WARN("Ensembl DB: " . $db_details->{'-DBNAME'});


    my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice');

    # my $exon = $ensembl_util->gene_adaptor->fetch_by_stable_id( $data->{ensembl_exon_id} );

    my $exon = $slice_adaptor->fetch_by_exon_stable_id(
        $data->{ensembl_exon_id}
    );

    my $exon_id = $exon->id;
    ## $exon_id

    my $exon_size = $exon->length;
    ## $exon_size
    my $exon_start = $exon->start;
    ## $exon_start
    my $exon_end = $exon->end;
    ## $exon_end


    $data->{exon_start} = $exon->start;
    $data->{exon_end} = $exon->end;


    ## $data
    find_target_crisprs( $data );

    undef $data;

# die;

}






=head2 find_target_crisprs

desc

=cut
sub find_target_crisprs {
    my ( $data ) = @_;

    my @crisprs;


    if ($data->{chr_strand} == 1 ) {
        $data->{insertion_site} = $data->{ins_start} + 3;
    } else {
        $data->{insertion_site} = $data->{ins_start} + 1;
    }




    @crisprs = crisprs_overlapping_region(
        {
            gene_id   => $data->{gene_id},
            chr_name  => $data->{chr_name},
            insertion_site => $data->{insertion_site},
            # chr_start => $data->{ins_start},
            # chr_end   => $data->{ins_end},
            strand    => $data->{chr_strand},
            # seq       => $data->{ins_seq},
        }
    );


    # unless ( @crisprs ) {
    #     $data->{fail_reason} = 'No crisprs found overlapping insertion site';
    #     # push @failed_targets, $data;
    #     # WARN( '.. no crisprs found' );
    #     return;
    # }

    # Check for valic (containing off-targets) crisprs, or let all crisprs be valid.
    # my @valid_crisprs = grep{ is_valid_crispr( $_ ) } map{ $_->as_hash } @crisprs;
    my @valid_crisprs = map{ $_->as_hash } @crisprs;

    unless ( @valid_crisprs ) {
        my $num_crisprs = @crisprs;
        $data->{fail_reason} = "Found $num_crisprs crispr(s) but non are valid";
        # push @failed_targets, $data;
        # WARN( ".. found $num_crisprs crisprs but non valid" );
        return;
    }
    # my $num_valid_crisprs = @valid_crisprs;
    # DEBUG( ".. found $num_valid_crisprs valid crisprs" );

    print_target( $data, \@valid_crisprs );

    return;
}



=head2 crisprs_overlapping_region

Diagram This...
Numbers allow insertion site to cut crispr site from between 2nd / 3rd bases and the 14th / 15th ones.

=cut
sub crisprs_overlapping_region {
    my ( $params ) = @_;

## $params

    # DEBUG("Getting crisprs for insertion site chr" . $params->{chr_name} . ':' . $params->{insertion_site} );
    print "Getting crisprs for gene " . $params->{gene_id} . " on insertion site chr" . $params->{chr_name} . ':' . $params->{insertion_site} . "\n";


my $insertion_site = $params->{insertion_site};
## $insertion_site

        ## REQUIRED >>> PAM out , last 4 bases out
        # my $pam_right_start = $params->{insertion_site} - 20;
        # my $pam_right_end   = $params->{insertion_site} - 4;
        # my $pam_left_start  = $params->{insertion_site} - 19;
        # my $pam_left_end    = $params->{insertion_site} - 3;

        my $pam_right_start = $params->{insertion_site} - 48;
        my $pam_right_end   = $params->{insertion_site} + 25;
        my $pam_left_start  = $params->{insertion_site} - 48;
        my $pam_left_end    = $params->{insertion_site} + 28;


## $pam_right_start


    my @pam_right_crisprs = $wge->resultset('Crispr')->search(
        {
            'species_id'  => $SPECIES_ID,
            'chr_name'    => $params->{chr_name},
            'pam_right'   => 1,
            'chr_start'   => {
                -between => [ $pam_right_start, $pam_right_end ],
            },
        },
    );

    my @pam_left_crisprs = $wge->resultset('Crispr')->search(
        {
            'species_id'  => $SPECIES_ID,
            'chr_name'    => $params->{chr_name},
            'pam_right'   => 0,
            'chr_start'   => {
                -between => [ $pam_left_start, $pam_left_end ],
            },
        },
    );

    return ( @pam_right_crisprs, @pam_left_crisprs );
}




=head2 print_target

Prints out the output file

=cut
sub print_target {
    my ( $data, $crisprs ) = @_;


    $data->{split_exon_L} = $data->{insertion_site} - $data->{exon_start} + 1;
    $data->{split_exon_R} = $data->{exon_end} - $data->{insertion_site};


    if ($data->{LAL_loc} && $data->{LAR_loc} && $data->{RAL_loc} && $data->{RAR_loc} ) {
        if ($data->{chr_strand} == 1 ) {

            $data->{LAL_loc} =~ /:(\d*)-/;
            my $la_start = $1;
            ## $la_start
            $data->{LAR_loc} =~ /-(\d*)$/;
            my $la_end = $1;
            ## $la_end
            $data->{RAL_loc} =~ /:(\d*)-/;
            my $ra_start = $1;
            ## $ra_start
            $data->{RAR_loc} =~ /-(\d*)$/;
            my $ra_end = $1;
            ## $ra_end

            $data->{LA_size} = $la_end - $la_start + 1;
            $data->{RA_size} = $ra_end - $ra_start + 1;



        } else {
            $data->{LAR_loc} =~ /:(\d*)-/;
            my $la_start = $1;
            ## $la_start
            $data->{LAL_loc} =~ /-(\d*)$/;
            my $la_end = $1;
            ## $la_end
            $data->{RAR_loc} =~ /:(\d*)-/;
            my $ra_start = $1;
            ## $ra_start
            $data->{RAL_loc} =~ /-(\d*)$/;
            my $ra_end = $1;
            ## $ra_end

            $data->{LA_size} = $la_end - $la_start + 1;
            $data->{RA_size} = $ra_end - $ra_start + 1;

        }
    } else {
        $data->{LA_size} = '';
        $data->{RA_size} = '';
    }


    my %base_params = (

        gene_id              => $data->{gene_id},
        marker_symbol        => $data->{marker_symbol},
        ensembl_gene_id      => $data->{ensembl_gene_id},
        ensembl_exon_id      => $data->{ensembl_exon_id},
        exon_size            => $data->{exon_size},
        exon_rank            => $data->{exon_rank},
        canonical_transcript => $data->{canonical_transcript},
        exon_start           => $data->{exon_start},
        exon_end             => $data->{exon_end},
        chr_name             => $data->{chr_name},
        chr_strand           => $data->{chr_strand},
        ins_start            => $data->{ins_start},
        ins_end              => $data->{ins_end},
        ins_seq              => $data->{ins_seq},
        insertion_site       => $data->{insertion_site},
        split_exon_L         => $data->{split_exon_L},
        split_exon_R         => $data->{split_exon_R},
        assembly             => $data->{assembly},
        build                => $data->{build},
        species              => $data->{species},
        LAL_seq              => $data->{LAL_seq},
        LAL_loc              => $data->{LAL_loc},
        LAL_gc               => $data->{LAL_gc},
        LAL_tm               => $data->{LAL_tm},
        LAR_seq              => $data->{LAR_seq},
        LAR_loc              => $data->{LAR_loc},
        LAR_gc               => $data->{LAR_gc},
        LAR_tm               => $data->{LAR_tm},
        LA_size              => $data->{LA_size},
        RAL_seq              => $data->{RAL_seq},
        RAL_loc              => $data->{RAL_loc},
        RAL_gc               => $data->{RAL_gc},
        RAL_tm               => $data->{RAL_tm},
        RAR_seq              => $data->{RAR_seq},
        RAR_loc              => $data->{RAR_loc},
        RAR_gc               => $data->{RAR_gc},
        RAR_tm               => $data->{RAR_tm},
        RA_size              => $data->{RA_size},

    );

    for my $crispr ( @{ $crisprs } ) {
        my %crispr_params = %base_params;

        my $proximity = $data->{insertion_site} - $crispr->{chr_start};
        ## $proximity

        $crispr_params{proximity}          = $proximity;
        $crispr_params{crispr_id}          = $crispr->{id};
        # $crispr_params{chr_name}         = $crispr->{chr_name};
        $crispr_params{crispr_start}       = $crispr->{chr_start};
        $crispr_params{crispr_end}         = $crispr->{chr_end};
        $crispr_params{crispr_seq}         = $crispr->{seq};
        $crispr_params{pam_right}          = $crispr->{pam_right};
        $crispr_params{off_target_summary} = $crispr->{off_target_summary};

        my ($off_zero, $off_one);
        if ($crispr->{off_target_summary}) {
            $crispr->{off_target_summary} =~ m/{0: (\d*), 1: (\d*), 2: (\d*), 3: (\d*), 4: (\d*)}/;
            $crispr_params{off_target_count} = $1 + $2 + $3 + $4 + $5;
            ($off_zero, $off_one) = ($1, $2);
        } else {
            $crispr_params{off_target_count} = '';
        }

## %crispr_params

        next if ($crispr_params{off_target_count} eq '');

        next if ($off_zero > 1 );

        # print "proximity: $proximity\n";


        if ($crispr->{pam_right}) {
            next if ( $proximity < 5 || $proximity > 22 );
        } else {
            next if ( $proximity < 1 || $proximity > 18 );
        }

        # print "print it!!!!\n";

        open($output_fh, ">>", "crisprs_${targets_file}");
        # print "[ @crispr_params{ @COLUMN_HEADERS } ]\n";
        $output_csv->print( $output_fh, [ @crispr_params{ @COLUMN_HEADERS } ] );
        close($output_fh);
        undef %crispr_params;


    }
    undef %base_params;

    return;
}