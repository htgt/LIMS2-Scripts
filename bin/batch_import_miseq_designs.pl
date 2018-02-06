#! /usr/bin/perl

use LIMS2::Model;
use strict;
use warnings;
use Try::Tiny;
use Carp;
use Getopt::Long;
use Data::Dumper;
use Moose;
use LIMS2::REST::Client;
use Bio::Perl qw( revcom );
use WGE::Model::DB;
use Text::CSV;
use WGE::Util::OffTargetServer;
use Path::Class;
use DesignCreate::Util::BWA;
use YAML::Any qw(DumpFile);
use JSON;
use POSIX qw(strftime);

use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);
my $logger = Log::Log4perl->get_logger('primer_design');

has lims2_api => (
    is         => 'ro',
    isa        => 'LIMS2::REST::Client',
    traits     => [ 'NoGetopt' ],
    lazy_build => 1
);



sub _build_lims2_api {
    my $self = shift;

    return LIMS2::REST::Client->new_with_config();
}

sub extract_rows {
    my ($csv, $fh, $data) = @_;

    my $headers = $csv->getline($fh);
    $csv->column_names( @{ $headers} );

    while (my $row = $csv->getline_hr($fh)) {
        unless($row->{WGE_ID} =~ /^\d+?$/) {
            print $row->{Experiment} . " - WGE Crispr ID is missing.\n";
            next;
        };
        $data->{$row->{Experiment}} = {
            wge_id      => $row->{WGE_ID},
            lab_crispr_id   => $row->{CRISPR_ID}, 
            crispr_seq  => $row->{'CRISPR Sequence'},
        };

        $data->{$row->{Experiment}}->{primers} = {
            exf => {
                seq     => $row->{'PCR forward'},
            },
            exr => {
                seq     => $row->{'PCR reverse'},
            },
            inf => {
                seq     => $row->{'MiSEQ forward'},
            },
            inr => {
                seq     => $row->{'MiSEQ reverse'},
            },
        };
    }

    return $data;
}

sub generate_bwa_query_file {
    my ($exp, $crispr, $data) = @_;

    my $root_dir = $ENV{ 'LIMS2_BWA_OLIGO_DIR' } // '/var/tmp/bwa';
    use Data::UUID;
    my $ug = Data::UUID->new();

    my $unique_string = $ug->create_str();
    my $dir_out = dir( $root_dir, '_' . $exp . '_' .  $unique_string );
    mkdir $dir_out->stringify  or die 'Could not create directory ' . $dir_out->stringify . ": $!";

    my $fasta_file_name = $dir_out->file( $exp . '_' . $crispr . '_oligos.fasta');
    my $fh = $fasta_file_name->openw();
    my $seq_out = Bio::SeqIO->new( -fh => $fh, -format => 'fasta' );

    foreach my $oligo ( sort keys %{ $data->{$exp}->{primers} } ) {
        my $fasta_seq = Bio::Seq->new( -seq => $data->{$exp}->{primers}->{$oligo}->{'seq'}, -id => $oligo );
        $seq_out->write_seq( $fasta_seq );
    }
    return ($fasta_file_name, $dir_out);
}

sub loci_builder {
    my ($oligo_hits, $exp, $oligo, $data, $strand) = @_;

    my $oligo_bwa = $oligo_hits->{$oligo};
    my $oligo_len = length($data->{$exp}->{primers}->{$oligo}->{seq});
    my $oligo_end = $oligo_bwa->{start} + $oligo_len;
    my $chr = $oligo_bwa->{chr};
    $chr =~ s/chr//;
    my $loci = {
        assembly    => 'GRCh38',
        chr_start   => $oligo_bwa->{start},
        chr_name    => $chr,
        chr_end     => $oligo_end,
        chr_strand  => $strand,
    };
    $data->{$exp}->{primers}->{$oligo}->{loci} = $loci;

    return $data;
}

sub yaml_oligos {
    my $yaml = shift;

    my $primers = $yaml->{primers};
    my @oligos;
    my $rev_oligo = {
        1   => {
            inf => 1,
            inr => -1,
            exf => 1,
            exr => -1,
        },
        -1  => {
            inf => -1,
            inr => 1,
            exf => -1,
            exr => 1,
        },
    };
    foreach my $primer (keys %$primers) {
        my $primer_data = $primers->{$primer};
        my $seq = $primer_data->{seq};
        if ($rev_oligo->{ $primer_data->{loci}->{chr_strand} }->{$primer} == -1) {
            $seq = revcom($seq)->seq;
        }
        my $oligo = {
            loci    => [ $primer_data->{loci} ],
            seq     => uc $seq,
            type    => uc $primer,
        };
        push(@oligos, $oligo);
    }

    $yaml->{oligos} = \@oligos;
    delete $yaml->{primers};

    return $yaml;
}

my ($file, $persist);

GetOptions(
    'file=s'  => \$file,
    'persist'   => \$persist,
)
or die usage_message();;

my $lims2_model = LIMS2::Model->new( user => 'lims2' );

my $date = strftime "%d-%m-%Y", localtime;
my $version = $lims2_model->software_version . '_' . $date;

my $design_parameters = {
    design_method       => 'miseq',
    'command-name'      => 'miseq-design-location',
    assembly            => 'GRCh38',
    created_by          => 'system',

    three_prime_exon    => 'null',
    five_prime_exon     => 'null',
    oligo_three_prime_align => '0',
    exon_check_flank_length =>  '0',
    primer_lowercase_masking    => 'null',
    num_genomic_hits            => "1",

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
    software_version    => $version,
};

my $data;
my $csv = Text::CSV->new({ binary => 1 }) or die "Cannot use CSV: " . Text::CSV->error_diag();

my $fh;
open $fh, "<:encoding(utf8)", $file or die;
$data = extract_rows($csv, $fh, $data);
close $fh;

my $ots_server = WGE::Util::OffTargetServer->new;

#Build data hash
foreach my $exp (keys %$data) {
    my @wge_crispr_arr = [$data->{$exp}->{wge_id}];
    my @crispr_arr = $lims2_model->import_wge_crisprs(\@wge_crispr_arr, 'Human', 'GRCh38');
    my $crispr_hash = $crispr_arr[0]->{db_crispr}->as_hash;
    $data->{$exp}->{lims_crispr} = $crispr_hash->{id};
    $data->{$exp}->{loci} = {
        assembly    => 'GRCh38',
        chr_start   => $crispr_hash->{locus}->{chr_start},
        chr_end     => $crispr_hash->{locus}->{chr_end},
        chr_name    => $crispr_hash->{locus}->{chr_name},
        chr_strand  => $crispr_hash->{locus}->{chr_strand},
    };
    my ($fasta, $dir) = generate_bwa_query_file($exp, $data->{$exp}->{wge_id}, $data);
    my $bwa = DesignCreate::Util::BWA->new(
            query_file        => $fasta,
            work_dir          => $dir,
            species           => 'Human',
            three_prime_check => 0,
            num_bwa_threads   => 2,
    );

    $bwa->generate_sam_file;
    my $oligo_hits = $bwa->oligo_hits;
    my $strand = 1;
    if ($oligo_hits->{exf}->{start} > $oligo_hits->{exr}->{start}) {
        $strand = -1;
    }
    foreach my $oligo (keys %$oligo_hits) {
        $data = loci_builder($oligo_hits, $exp, $oligo, $data, $strand);
    }
}

my $json = JSON->new->allow_nonref;
my $base = $ENV{YAML_DUMP};
#YAML creation
foreach my $experiment (keys %$data) {
    my $design_yaml = $data->{$experiment};
    my $params = $design_parameters;

    my @lab_id = split(/\_/, $data->{$experiment}->{lab_crispr_id});
    my $symbol = $lab_id[0];

    my $search_terms = {
        species     => 'Human',
        search_term => $symbol,
    };
    my $gene = $lims2_model->find_gene($search_terms);

    $params->{species} = 'Human';
    $params->{chr_name} = $gene->{chromosome};
    $params->{chr_strand} = $design_yaml->{loci}->{chr_strand};
    $params->{target_start} = $design_yaml->{loci}->{chr_start};
    $params->{target_end} = $design_yaml->{loci}->{chr_end};
    $params->{target_genes} = [ $gene->{gene_id} ];

    my $params_json = $json->encode($design_parameters);
    $design_yaml->{design_parameters} = $params_json;
    $design_yaml->{created_by} = 'system';
    $design_yaml->{species} = 'Human';
    $design_yaml->{type} = 'miseq';
    $design_yaml->{name} = $experiment;
    $design_yaml->{gene_ids} = [{
        gene_id => $gene->{gene_id},
        gene_type_id => 'HGNC',
    }];

    $design_yaml = yaml_oligos($design_yaml);
    my $lims_crispr = $design_yaml->{lims_crispr};

    delete $design_yaml->{crispr_seq};
    delete $design_yaml->{lab_crispr_id};
    delete $design_yaml->{loci};
    delete $design_yaml->{wge_id};
    delete $design_yaml->{lims_crispr};
    
    if ($base) {
        my $file_name = $base . '/' . $experiment . '-' . $lims_crispr . '.yaml';
        my $fh;
        open ($fh, '>', $file_name) or die "$!";
        DumpFile($file_name, $design_yaml);
        close $fh;
    }
}

if ($persist) {
    my $lims = {
        lims2_api         => _build_lims2_api(),
        design_method     => 'miseq',
    };
    my $self = $lims->{lims2_api}->_build_ua();
    
    opendir my $dirh, $base or die "Can't open dir $base: $!";
    my @yamls = sort grep {/\.yaml/} readdir $dirh;
    closedir $dirh;

    my $csv = Text::CSV->new({binary => 1, eol => "\n"}) or die "Cannot use CSV: ".Text::CSV->error_diag ();
    my $res_file_path = $base . '/results.csv';

    open my $res, '>', $res_file_path or die "$res_file_path: $!";
    my $header = ['Experiment','Design ID','LIMS2 Crispr'];
    $csv->print($res, $header);
    foreach my $yaml (@yamls) {
        my $design_data = YAML::Any::LoadFile( $base . '/' . $yaml );
        my $design = $lims->{lims2_api}->POST('design', $design_data );
        print ('Design persisted: ' . $design->{id} . ", Experiment YAML file: " . $yaml  ."\n");

        my $file_name = $yaml;
        $file_name =~ s/.yaml//;
        my @details = split /\-/, $file_name;
        my $row = [$details[0], $design->{id}, $details[1]];
        $csv->print($res, $row);
    }
    close $res;
}

1;
