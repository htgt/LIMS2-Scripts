#! /usr/bin/perl

use LIMS2::Model;
use feature qw/ say /;
use strict;
use warnings;
use Try::Tiny;


use Bio::EnsEMBL::Registry;
use Smart::Comments;

use lib '/nfs/users/nfs_t/tg6/tg6misc/Design-Creation/lib';
use DesignCreate::Util::Primer3;
use LIMS2::Model::Util::OligoSelection;

use Bio::SeqIO;
use Path::Class;

my $model = LIMS2::Model->new( { user => 'webapp', audit_user => $ENV{USER} .'@sanger.ac.uk' } );

print "starting...\n";

# csv file with a list of targets (inserion sites)
my $targets_file = shift @ARGV;

# the format of the file is as above (header on first row)
# gene_id marker_symbol ensembl_gene_id ensembl_exon_id exon_size exon_rank canonical_transcript chr_name chr_strand ins_start ins_end  ins_seq assembly build species
# HGNC:5  A1BG          ENSG00000121410 ENSE00003598943 270       3         ENST00000263100      19       -1         58353117  58353120 CAGG    GRCh38   79    Human







# get ensembl ready to query
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous');





open ( my $fh, '<', $targets_file ) or die( "Can not open $targets_file " . $! );

open ( OUTPUT, '>', "primers_$targets_file" ) or die( "Can not open output file " . $! );


my $header = <$fh>;
chomp $header;
# ### $header
# my @header_fields = split(',', $header);

# ### @header_fields

print OUTPUT "$header,LAL_seq,LAL_loc,LAL_gc,LAL_tm,LAR_seq,LAR_loc,LAR_gc,LAR_tm,RAL_seq,RAL_loc,RAL_gc,RAL_tm,RAR_seq,RAR_loc,RAR_gc,RAR_tm,\n";

while(<$fh>) {
    my $line = $_;
    chomp $line;
    my ( $gene_id, $marker_symbol, $ensembl_gene_id, $ensembl_exon_id, $exon_size, $exon_rank, $canonical_transcript, $chr_name, $chr_strand, $ins_start, $ins_end, $ins_seq, $assembly, $build, $species) =  split(',', $line);

    ## $line

    ## $ins_seq

    ## $chr_name
    ## $chr_strand
    ## $ins_start
    ## $ins_end

    # In this case, the target sequence is always in the format NAGX (N can be any base, X can be A or G). The cut site is between thr AG and the A or G. This this means an offset of 3 bases from the start position of the sequence.
    # $ins_seq =~ /(^.?AG)(A|G)/;
    # my $ga_site = $1;
    # my $offset = length($ga_site);


    my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice');



    # $insert_coord is the coordinate number of the base just after insertion site...
    my $insert_coord = $ins_start + 3;

    if ( $chr_strand == -1 ) {
        $insert_coord = $ins_start + 1;
    }

    ## $insert_coord



    # my $test_slice_region = $slice_adaptor->fetch_by_region(
    #     'chromosome',
    #     $chr_name,
    #     $insert_coord - 10,
    #     $insert_coord + 10,
    #     $chr_strand,
    # );

    # my $test_sequence = $test_slice_region->seq;
    # ### $test_sequence

    # my $test_chr_region_start = $test_slice_region->start;
    # ### $test_chr_region_start

    # my $test_chr_region_end = $test_slice_region->end;
    # ### $test_chr_region_end








    my $ra_start_coord =  $ins_start;
    my $ra_end_coord =  $ins_start + 1200 + 3;

    my $la_start_coord =  $ins_end - 1200 - 1;
    my $la_end_coord =  $ins_end;

    if ( $chr_strand == -1 ) {

        $ra_start_coord = $ins_end - 1200 - 3;
        $ra_end_coord =  $ins_end;

        $la_start_coord = $ins_start;
        $la_end_coord =  $ins_start + 1200 + 1;
    }

## $ra_start_coord
## $ra_end_coord

## $la_start_coord
## $la_end_coord

    my $ra_slice_region = $slice_adaptor->fetch_by_region(
        'chromosome',
        $chr_name,
        $ra_start_coord,
        $ra_end_coord,
        $chr_strand,
    );

    my $ra_sequence = $ra_slice_region->seq;
    ## $ra_sequence





    # my $la_start_coord =  $ins_start + $offset;
    # my $la_end_coord =  $la_start_coord + 1200;
    my $la_slice_region = $slice_adaptor->fetch_by_region(
        'chromosome',
        $chr_name,
        $la_start_coord,
        $la_end_coord,
        $chr_strand,
    );

    my $la_sequence = $la_slice_region->seq;
    ## $la_sequence


## RIGHT ARM

my $ra_region_bio_seq = Bio::Seq->new( -alphabet => 'dna', -seq => $ra_slice_region->seq, -verbose => -1 );
## $ra_region_bio_seq

    my $ra_p3 = DesignCreate::Util::Primer3->new_with_config(
        configfile => 'primer3_crispr_conditional_ra_config.yaml',

    );
    my $ra_dir_out = dir( '/nfs/users/nfs_t/tg6/tg6misc/LIMS2-Scripts' );
    my $ra_logfile = $ra_dir_out->file( 'primers.log');


    my ( $ra_result, $ra_primer3_explain ) = $ra_p3->run_primer3( $ra_logfile->absolute, $ra_region_bio_seq, # bio::seqI
            {
                # SEQUENCE_TARGET => '30,800',
            } );

### $ra_result

### $ra_primer3_explain

my $ra_primer_data;

    if ( $ra_result->num_primer_pairs ) {
        # INFO ( "primer pairs: " . $ra_result->num_primer_pairs );
        $ra_primer_data = parse_primer3_results( $ra_result );
        $ra_primer_data->{'error_flag'} = 'pass';
        # $primer_passes = genomic_check( $design_id, $well_id, $species, $primer_data, $chr_strand );
        # $primer_passes->{'genomic_error_flag'} = $primer_passes->{'pair_count'} > 0 ? 'pass' : 'fail';
    }

## $ra_primer_data

my ($ral_seq, $ral_loc, $ral_gc, $ral_tm, $rar_seq, $rar_loc, $rar_gc, $rar_tm ) = ('','','','','','','','');
my ($ra_left_start, $ra_left_end, $ra_right_start, $ra_right_end);


if ($ra_primer_data) {

$ra_left_start = $ra_start_coord + $ra_primer_data->{'left'}->{'left_0'}->{'location'}->start - 1;
$ra_left_end = $ra_start_coord + $ra_primer_data->{'left'}->{'left_0'}->{'location'}->end - 1;

$ra_right_start = $ra_start_coord + $ra_primer_data->{'right'}->{'right_0'}->{'location'}->start - 1;
$ra_right_end = $ra_start_coord + $ra_primer_data->{'right'}->{'right_0'}->{'location'}->end - 1;

# my $leftstart = $ra_primer_data->{'left'}->{'left_0'}->{'location'}->start;
# ### $leftstart
# my $leftend = $ra_primer_data->{'left'}->{'left_0'}->{'location'}->end;
# ### $leftend

# my $rightstart = $ra_primer_data->{'right'}->{'right_0'}->{'location'}->start;
# ### $rightstart
# my $rightend = $ra_primer_data->{'right'}->{'right_0'}->{'location'}->end;
# ### $rightend

    if ( $chr_strand == -1 ) {

        $ra_left_start = $ra_end_coord + 1 - $ra_primer_data->{'left'}->{'left_0'}->{'location'}->end;
        $ra_left_end = $ra_end_coord + 1 - $ra_primer_data->{'left'}->{'left_0'}->{'location'}->start;
        # $ra_left_start = $ra_start_coord - $ra_primer_data->{'left'}->{'left_0'}->{'location'}->end + $ra_primer_data->{'left'}->{'left_0'}->{'location'}->start;
        # $ra_left_end = $ra_start_coord;

        $ra_right_start = $ra_end_coord + 1 - $ra_primer_data->{'right'}->{'right_0'}->{'location'}->end;
        $ra_right_end = $ra_end_coord + 1 - $ra_primer_data->{'right'}->{'right_0'}->{'location'}->start;
        # $ra_right_start = $ra_end_coord - $ra_primer_data->{'right'}->{'right_0'}->{'location'}->end + $ra_primer_data->{'right'}->{'right_0'}->{'location'}->start;
        # $ra_right_end = $ra_end_coord;

    }


    $ral_seq = $ra_primer_data->{'left'}->{'left_0'}->{'seq'};
    $ral_loc = 'chr' . $chr_name . ':' . $ra_left_start . '-' . $ra_left_end;
    $ral_gc = $ra_primer_data->{'left'}->{'left_0'}->{'gc_content'};
    $ral_tm = $ra_primer_data->{'left'}->{'left_0'}->{'melting_temp'};
    $rar_seq = $ra_primer_data->{'right'}->{'right_0'}->{'seq'};
    $rar_loc = 'chr' . $chr_name . ':' . $ra_right_start . '-' . $ra_right_end;
    $rar_gc = $ra_primer_data->{'right'}->{'right_0'}->{'gc_content'};
    $rar_tm = $ra_primer_data->{'right'}->{'right_0'}->{'melting_temp'};

}


## $rar_seq
## $rar_loc

## $ral_seq
## $ral_loc










## LEFT ARM

my $la_region_bio_seq = Bio::Seq->new( -alphabet => 'dna', -seq => $la_slice_region->seq, -verbose => -1 );
## $la_region_bio_seq

    my $la_p3 = DesignCreate::Util::Primer3->new_with_config(
        configfile => 'primer3_crispr_conditional_la_config.yaml',

    );
    my $la_dir_out = dir( '/nfs/users/nfs_t/tg6/tg6misc/LIMS2-Scripts' );
    my $la_logfile = $la_dir_out->file( 'primers.log');


    my ( $la_result, $la_primer3_explain ) = $la_p3->run_primer3( $la_logfile->absolute, $la_region_bio_seq, # bio::seqI
            {
                # SEQUENCE_TARGET => '30,800',
            } );

### $la_result

### $la_primer3_explain

my $la_primer_data;

    if ( $la_result->num_primer_pairs ) {
        # INFO ( "primer pairs: " . $la_result->num_primer_pairs );
        $la_primer_data = parse_primer3_results( $la_result );
        $la_primer_data->{'error_flag'} = 'pass';
        # $primer_passes = genomic_check( $design_id, $well_id, $species, $primer_data, $chr_strand );
        # $primer_passes->{'genomic_error_flag'} = $primer_passes->{'pair_count'} > 0 ? 'pass' : 'fail';
    }

## $la_primer_data


my ($lal_seq, $lal_loc, $lal_gc, $lal_tm, $lar_seq, $lar_loc, $lar_gc, $lar_tm ) = ('','','','','','','','');
my ($la_left_start, $la_left_end, $la_right_start, $la_right_end);

if ($la_primer_data) {


$la_left_start = $la_start_coord - 1 + $la_primer_data->{'left'}->{'left_0'}->{'location'}->start;
$la_left_end = $la_start_coord - 1 + $la_primer_data->{'left'}->{'left_0'}->{'location'}->end;

$la_right_start = $la_start_coord - 1 + $la_primer_data->{'right'}->{'right_0'}->{'location'}->start;
$la_right_end = $la_start_coord - 1 + $la_primer_data->{'right'}->{'right_0'}->{'location'}->end;

# my $leftstart = $la_primer_data->{'left'}->{'left_0'}->{'location'}->start;
# ### $leftstart
# my $leftend = $la_primer_data->{'left'}->{'left_0'}->{'location'}->end;
# ### $leftend

# my $rightstart = $la_primer_data->{'right'}->{'right_0'}->{'location'}->start;
# ### $rightstart
# my $rightend = $la_primer_data->{'right'}->{'right_0'}->{'location'}->end;
# ### $rightend


    if ( $chr_strand == -1 ) {

        $la_left_start = $la_end_coord + 1 - $la_primer_data->{'left'}->{'left_0'}->{'location'}->end;
        $la_left_end = $la_end_coord + 1 - $la_primer_data->{'left'}->{'left_0'}->{'location'}->start;

        $la_right_start = $la_end_coord + 1 - $la_primer_data->{'right'}->{'right_0'}->{'location'}->end;
        $la_right_end = $la_end_coord + 1 - $la_primer_data->{'right'}->{'right_0'}->{'location'}->start;

    }

    $lal_seq = $la_primer_data->{'left'}->{'left_0'}->{'seq'};
    $lal_loc = 'chr' . $chr_name . ':' . $la_left_start . '-' . $la_left_end;
    $lal_gc = $la_primer_data->{'left'}->{'left_0'}->{'gc_content'};
    $lal_tm = $la_primer_data->{'left'}->{'left_0'}->{'melting_temp'};
    $lar_seq = $la_primer_data->{'right'}->{'right_0'}->{'seq'};
    $lar_loc = 'chr' . $chr_name . ':' . $la_right_start . '-' . $la_right_end;
    $lar_gc = $la_primer_data->{'right'}->{'right_0'}->{'gc_content'};
    $lar_tm = $la_primer_data->{'right'}->{'right_0'}->{'melting_temp'};
}
## $lal_seq
## $lal_loc
## $lar_seq
## $lar_loc
















# die;


print OUTPUT "$line,$lal_seq,$lal_loc,$lal_gc,$lal_tm,$lar_seq,$lar_loc,$lar_gc,$lar_tm,$ral_seq,$ral_loc,$ral_gc,$ral_tm,$rar_seq,$rar_loc,$rar_gc,$rar_tm\n";

print "$line,$lal_seq,$lal_loc,$lal_gc,$lal_tm,$lar_seq,$lar_loc,$lar_gc,$lar_tm,$ral_seq,$ral_loc,$ral_gc,$ral_tm,$rar_seq,$rar_loc,$rar_gc,$rar_tm\n";



}

















sub get_repeat_masked_sequence {
    my $params = shift;

    my $slice_region = $params->{'slice_region'};
    my $repeat_mask = $params->{'repeat_mask'};
    my $revcom = $params->{'revcom'};
    my $seq;
    if ( $repeat_mask->[0] eq 'NONE' ) {
        DEBUG('No repeat masking selected');
        $seq = Bio::Seq->new( -alphabet => 'dna', -seq => $slice_region->seq, -verbose => -1 );
    }
    else {
        DEBUG('Repeat masking selected');
        $seq = Bio::Seq->new( -alphabet => 'dna', -seq => $slice_region->get_repeatmasked_seq($repeat_mask)->seq, -verbose => -1 );
    }
    if ( $revcom ) {
        $seq = $seq->revcom;
    }
    return $seq;
}


sub parse_primer3_results {
    my $result = shift;

    my $oligo_data;
    # iterate through each primer pair
    $oligo_data->{pair_count} = $result->num_primer_pairs;
    while (my $pair = $result->next_primer_pair) {
        # do stuff with primer pairs...
        my ($fp, $rp) = ($pair->forward_primer, $pair->reverse_primer);
        $oligo_data->{'left'}->{$fp->display_name} = parse_primer( $fp );
        $oligo_data->{'right'}->{$rp->display_name} = parse_primer( $rp );
    }

    return $oligo_data;
}

sub parse_primer {
    my $primer = shift;

    my %oligo_data;

    my @primer_attrs = qw/
        length
        melting_temp
        gc_content
        rank
        location
    /;


    %oligo_data = map { $_  => $primer->$_ } @primer_attrs;
    $oligo_data{'seq'} = $primer->seq->seq;

    return \%oligo_data;
}