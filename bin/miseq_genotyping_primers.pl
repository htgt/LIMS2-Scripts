#! /usr/bin/perl

use LIMS2::Model;
use strict;
use warnings;
use Try::Tiny;
use Carp;
use Getopt::Long;
use Data::Dumper;
use Moose;
use Bio::Perl qw( revcom );
use WGE::Model::DB;
use LIMS2::Model::Util::OligoSelection qw(
        pick_crispr_primers
        pick_single_crispr_primers
        pick_miseq_internal_crispr_primers
        pick_miseq_crispr_PCR_primers
        oligo_for_single_crispr
        pick_crispr_PCR_primers
);

use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);
my $logger = Log::Log4perl->get_logger('primer_design');

has lims2_api => (
    is         => 'ro',
    isa        => 'LIMS2::REST::Client',
    traits     => [ 'NoGetopt' ],
    lazy_build => 1
);

my ($crispr, $species, $wge, $persist);

GetOptions(
    'crispr=s'  => \$crispr,
    'species=s' => \$species,
    'wge'       => \$wge,
    'persist'   => \$persist,
)
or die usage_message();;

sub usage_message {
return << "END_DIE";
Usage: perl Miseq_genotyping_primers.pl
    --crispr = LIMS2 / WGE crispr ID
    --species = Human
    --wge = Toggle if using a WGE crispr
END_DIE
}

sub _build_lims2_api {
    my $self = shift;

    return LIMS2::REST::Client->new_with_config();
}

sub find_appropriate_primers {
    my ($crispr_primers, $target, $max, $region, $crispr) = @_;

    my @primers = keys %{$crispr_primers->{left}};
    my $closest->{record} = 5000;
    my @test;
    foreach my $prime (@primers) {
        my $int = (split /_/, $prime)[1];
        my $left_location_details = $crispr_primers->{left}->{'left_' . $int}->{location};
        my $right_location_details = $crispr_primers->{right}->{'right_' . $int}->{location};
        my $range = $right_location_details->{_start} - $left_location_details->{_end};
        my $start_coord = index ($region, $crispr_primers->{left}->{'left_' . $int}->{seq});
        my $end_coord = index ($region, revcom($crispr_primers->{right}->{'right_' . $int}->{seq})->seq);
        my $primer_diff = abs (($end_coord - 1022) - (1000 - $start_coord));

        my $primer_range = {
            name    => '_' . $int,
            start   => $left_location_details->{_end},
            end     => $right_location_details->{_start},
            lseq    => $crispr_primers->{left}->{'left_' . $int}->{seq},
            rseq    => $crispr_primers->{right}->{'right_' . $int}->{seq},
            range   => $range,
            diff    => $primer_diff,
        };

        push @test, $primer_range;
        if ($range < $max) {
            my $amplicon_score = ($target - $range) + $primer_diff;
            if ($amplicon_score < $closest->{record}) {
                $closest = {
                    record  => $amplicon_score,
                    primer  => $int,
                };
            }
        }
    }
    print Dumper $closest;
    return $crispr_primers->{left}->{'left_' . $closest->{primer}}, $crispr_primers->{right}->{'right_' . $closest->{primer}};
}
my $wge_model = WGE::Model::DB->new;
my $lims2_model = LIMS2::Model->new( user => 'lims2' );

my $crispr_id;
if ($wge) {
    my @wge_crispr_arr;
    push (@wge_crispr_arr, $crispr);
    my @crispr_arr = $lims2_model->import_wge_crisprs(\@wge_crispr_arr, $species, 'GRCh38');
    $crispr_id = $crispr_arr[0]->{lims2_id};
} else {
    $crispr_id = $crispr;
}
my $params = {
    crispr_id => $crispr_id,
    species => $species,
    repeat_mask => [''],
    offset => 20,
    increment => 15,
    well_id => 'Miseq_Crispr_' . $crispr,
};

$ENV{'LIMS2_SEQ_SEARCH_FIELD'} = "170";
$ENV{'LIMS2_SEQ_DEAD_FIELD'} = "50";

my ($internal_crispr, $internal_crispr_primers) = pick_miseq_internal_crispr_primers($lims2_model, $params);

$ENV{'LIMS2_PCR_SEARCH_FIELD'} = "350";
$ENV{'LIMS2_PCR_DEAD_FIELD'} = "170";
$params->{increment} = 50;

my $crispr_seq = {
    chr_region_start    => $internal_crispr->{left_crispr}->{chr_start},
    left_crispr         => { chr_name   => $internal_crispr->{left_crispr}->{chr_name} },
};
my $en_strand = {
    1   => 'plus',
    -1  => 'minus',
};
$DB::single=1;
if ($internal_crispr_primers->{error_flag} eq 'fail' ) {
    print "Primer generation failed: Internal primers - " . $internal_crispr_primers->{error_flag} . "; Crispr:" . $crispr . "\n";
    exit;
}
$params->{crispr_primers} = { 
    crispr_primers  => $internal_crispr_primers,
    crispr_seq      => $crispr_seq,
    strand          => $en_strand->{$internal_crispr->{left_crispr}->{chr_strand}},
};

my ($pcr_crispr, $pcr_crispr_primers) = pick_miseq_crispr_PCR_primers($lims2_model, $params);

$DB::single=1;
if ($pcr_crispr_primers->{error_flag} eq 'fail') {
    print "Primer generation failed: PCR results - " . $pcr_crispr_primers->{error_flag} . "; Crispr " . $crispr . "\n";
    exit;
} elsif ($pcr_crispr_primers->{genomic_error_flag} eq 'fail') {
    print "PCR genomic check failed; PCR results - " . $pcr_crispr_primers->{genomic_error_flag} . "; Crispr " . $crispr . "\n";
    exit;
}

my $slice_adaptor = $lims2_model->ensembl_slice_adaptor($species);
my $slice_region = $slice_adaptor->fetch_by_region(
    'chromosome',
    $internal_crispr->{'left_crispr'}->{'chr_name'},
    $internal_crispr->{'left_crispr'}->{'chr_start'} - 1000,
    $internal_crispr->{'left_crispr'}->{'chr_end'} + 1000,
    1,
);
my $crispr_loc = index ($slice_region->seq, $internal_crispr->{left_crispr}->{seq});
my ($inf, $inr) = find_appropriate_primers($internal_crispr_primers, 260, 297, $slice_region->seq, $crispr_loc);
my ($exf, $exr) = find_appropriate_primers($pcr_crispr_primers, 750, 3000, $slice_region->seq);

#BWA

my $primers = {
    crispr  => $internal_crispr->{left_crispr}->{seq},
    inf     => $inf->{seq},
    inr     => $inr->{seq},
    exf     => $exf->{seq},
    exr     => $exr->{seq},
};
my $order = {
    crispr  => 3,
    inf     => 2,
    inr     => 4,
    exf     => 1,
    exr     => 5,
};

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
foreach my $key (keys %$primers) {
    if ($key eq 'crispr') {
        print ">3-$crispr_id-$key\n";
    } else {
        print ">$order->{$key}-$crispr_id-$key\n";
    }
    print "$primers->{$key}\n";
}

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
my $loc = $internal_crispr->{left_crispr};
if ($persist) {
    my $design_parameters = {
        design_method       => 'miseq',
        'command-name'      => 'miseq-design-location',
        assembly            => 'GRCh38',
        species             => $species,
        created_by          => 'system',
        chr_name            => $loc->{chr_name},
        chr_strand          => $loc->{chr_strand},
        target_genes        => ['X'],

        three_prime_exon    => 'null',
        five_prime_exon     => 'null',
        oligo_three_prime_align => '0',
        exon_check_flank_length =>  '0',
        primer_lowercase_masking    => 'null',
        num_genomic_hits            => "1",

        target_start        => $loc->{chr_start},
        target_end          => $loc->{chr_end},

        region_length_3F    => '20',
        region_length_3R    => '20',
        region_length_5F    => '20',
        region_length_5R    => '20',

        region_offset_3F    => 80,
        region_offset_3R    => 80,
        region_offset_5F    => 300,
        region_offset_5R    => 300,

        primer_min_size     => '18',
        primer_min_gc       => '40',
        primer_opt_gc_content   => '50',
        primer_opt_size     => '20',
        primer_max_size     => '22',
        primer_max_gc       => '60',
       
        primer_min_tm       => '57',
        primer_opt_tm       => '60',
        primer_max_tm       => '63',

        repeat_mask_class   => [],
        
        'ensembl-version'   => 'X',
        software_version    => 'X',
    };

    my $design_data = {
        created_by  => 'system',
        design_parameters   => $design_parameters,
    };
    my $dir; 
    my $lims = {
        lims2_api         => _build_lims2_api(),
        dir               => $dir,
        design_method     => 'conditional-inversion',
    };

    my $self = $lims->{lims2_api}->_build_ua();
    my $design = $lims->{lims2_api}->POST('design', {'z' => 'y'} );
    print ('Design persisted: ' . $design->{id} );
}

1;
