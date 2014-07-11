#! /usr/bin/perl

use LIMS2::Model;
use strict;
use warnings;
use Try::Tiny;
use Carp;
use Getopt::Long;

use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);
my $logger = Log::Log4perl->get_logger('primer_design');

use LIMS2::Model::Util::OligoSelection qw(
        pick_crispr_primers
        pick_single_crispr_primers
);

my $plate_name_param = '';
my $plate_well_param = '';
my $species_param = 'Human';
my $assembly_param = '';
my $crispr_type = 'pair';
my $format_param = '';
my @repeat_mask_param;
my @d_plate_param;
my $crispr_pair_id_param = '';
my $genotyping_primers_required = 'YES';
my $crispr_primers_required = 'YES';
my $pcr_primers_required = 'YES';
my $persist_param = 'BOTH';
my $plate_crispr_left_param;
my $plate_crispr_right_param;

GetOptions(
    'plate=s'       => \$plate_name_param,
    'well=s'        => \$plate_well_param,
    'repeat_mask=s' => \@repeat_mask_param,
    'species=s'     => \$species_param,
    'crispr_type=s' => \$crispr_type,
    'format=s'      => \$format_param,
    'assembly=s'    => \$assembly_param,
    'd_plate=s'     => \@d_plate_param,
    'pair_id=s'     => \$crispr_pair_id_param,
    'genotyping_primers=s'  => \$genotyping_primers_required,
    'crispr_primers=s'  => \$crispr_primers_required,
    'pcr_primers=s'  => \$pcr_primers_required,
    'persist=s'     => \$persist_param,
    'crispr_left=s'   => \$plate_crispr_left_param,
    'crispr_right=s'  => \$plate_crispr_right_param,
)
or die usage_message();;

if ($plate_name_param eq '') {
    die usage_message()
}

if ( scalar @repeat_mask_param == 0 ) {
    push  @repeat_mask_param,'NONE';
}

$persist_param = uc($persist_param);

my $lims2_model = LIMS2::Model->new( { user => 'lims2' } );
my $lims2_schema = $lims2_model->schema;
if ($assembly_param eq '') {
    my $assembly_r = $lims2_schema->resultset('SpeciesDefaultAssembly')->find( { species_id => $species_param } );
    $assembly_param = $assembly_r->assembly_id;
}

INFO( '-------------------- New Primer Calculation Begins -----------------------');
my $data_clip;
if ( $plate_well_param  eq '') {
    if ( $crispr_type eq 'single' ) {
        $data_clip = primers_for_single_crispr_plate({
            'model' => $lims2_model,
            'plate_name' => $plate_name_param,
            'repeat_mask_list' => \@repeat_mask_param,
            'd_plate_list' => \@d_plate_param,
            'species' => $species_param,
            'assembly' => $assembly_param,
        });
    }
    else {
        $data_clip = primers_for_plate({
            'model' => $lims2_model,
            'plate_name' => $plate_name_param,
            'plate_crispr_right' => $plate_crispr_right_param,
            'plate_crispr_left' => $plate_crispr_left_param,
            'repeat_mask_list' => \@repeat_mask_param,
            'species' => $species_param,
            'assembly' => $assembly_param,
        });
    }
}
else
{
    if ( $crispr_type eq 'single' ) {
        # $plate_well_param is optional here
        $data_clip = primers_for_single_crispr_plate({
            'model' => $lims2_model,
            'plate_name' => $plate_name_param,
            'repeat_mask' =>\@repeat_mask_param,
            'd_plate_list' => \@d_plate_param,
            'species' => $species_param,
            'assembly' => $assembly_param,
            'plate_well' => $plate_well_param,
        });
    }
    else {
        $data_clip = primers_for_plate({
            'model' => $lims2_model,
            'plate_name' => $plate_name_param,
            'plate_crispr_right' => $plate_crispr_right_param,
            'plate_crispr_left' => $plate_crispr_left_param,
            'repeat_mask_list' => \@repeat_mask_param,
            'species' => $species_param,
            'assembly' => $assembly_param,
            'plate_well' => $plate_well_param,
            'crispr_pair_id' => $crispr_pair_id_param,
        });
    }
}

# If the format is specified as fsa, fasta format files will be generated for each well suitable for input to BlastN
if ( $format_param eq 'fsa' ) {
# Checkout $data_clip
    $logger->info('Generating blastable fasta files');
    my $fsa_rows = generate_fsa_files( $data_clip );
}

exit;

sub usage_message {
return << "END_DIE";

Usage: perl genotyping_primers.pl
    --plate=plate_name
    [--crispr_right=crispr_right_plate_name]
    [--crispr_left=crispr_left_plate_name]
    [--well=well_name]
    [--crispr_type=(single | pair)]
    [--species=(Human | Mouse)]
    [--assembly=assembly_name (e.g., GRCm38), default is the LIMS2 default assembly
    [--repeat_mask=TRF [--repeat_mask=...] (default: NONE)]
    [--format=fsa]
    [--d_plate=plate_name [--d_plate=...]] (default: no plates) These plates hold a list of designs used for multi-design disambiguation
    [--genotyping_primers=[YES | NO] (default: YES)
    [--crispr_primers=[YES | NO] (default: YES)
    [--pcr_primers=[YES | NO] (default: YES)
    [--persist= [file | db | both] (default: both)

Optional parameters in square brackets
Default species is Human
Default crispr_type is pair
If you do not supply crispr_left/crispr_right plates, the script tries to ascend the crispr_v hierarchy to find the crispr_pair_ids

END_DIE
}

# subroutines start here
#
sub primers_for_single_crispr_plate {
    my $p = shift;

    my $model = $p->{'model'};
    my $plate_name_input = $p->{'plate_name'};
    my $repeat_mask = $p->{'repeat_mask'};
    my $d_plate_list = $p->{'d_plate_list'};
    my $species_input = $p->{'species'};
    my $assembly_input = $p->{'assembly'};
    my $plate_well_input = $p->{'plate_well'};

    my $dis_designs; # hash ref for disambiguation design list

    $logger->info("Starting primer generation for plate $plate_name_input");
    $logger->info("Primers for sequencing and pcr of a single crispr per well" );
    $logger->info("Using species: $species_input (assembly: $assembly_input)" );
    my $rpt_string = join( ',', @$repeat_mask);
    $logger->info("Using sequence repeat mask of: $rpt_string");
    if ( scalar @$d_plate_list ) {
        $logger->info('Disambiguation plate list defined:');
        my $counter = 0;
        foreach my $d_plate_name ( @$d_plate_list ){
            $counter++;
            $logger->info($counter . ': ' . $d_plate_name);
        }
        $dis_designs = create_d_plate_hash( $model, $d_plate_list );
    }


    my $species = $species_input;
    my $assembly = $assembly_input;

    my $plate_rs = $model->schema->resultset( 'Plate' )->search(
        {
            'name' => $plate_name_input
        },
    );

     if ( !$plate_rs->count ) {
         $logger->fatal("No wells on plate or non-existent plate: $plate_name_input");
         exit;
     }


    my $plate = $plate_rs->first;
    my $plate_name = $plate->name;
    $logger->info( 'Plate selected: ' . $plate_name );

    my @wells = $plate->wells->all;

    my $well_count = @wells;
    $logger->info( 'Processing crispr primers for ' . $well_count . ' wells:');
    my @well_id_list;

    # Select just the well we want
    if ( $plate_well_input ) {
        $logger->info("Selecting well: [$plate_well_input]");
        @wells = grep { $_->name =~ /$plate_well_input/  } @wells;
    }

    foreach my $well ( @wells ) {
        push @well_id_list, $well->id;
    }

    my $design_data_cache = $model->create_design_data_cache(
            \@well_id_list,
        );

    if (! $design_data_cache ) {
        $logger->debug( 'Unable to populate design data cache (no designs) - deferring to crispr phase' );
    }
    else
    {
        $logger->debug( 'Created design well data cache' );
    }

    my $sc_params = {
            'wells' => \@wells,
            'design_data_cache' => $design_data_cache,
            'model' => $model,
            'species' => $species,
            'assembly' => $assembly,
            'plate_name' => $plate_name,
            'repeat_mask' => $repeat_mask,
            'dis_designs' => $dis_designs,
        };

    if ( $plate_well_input ) {
        $sc_params->{'plate_well'} = $plate_well_input;
    }

    my $clip = run_single_crispr_primers( $sc_params );

    $logger->info( 'Done' );
    return $clip;
}

sub primers_for_plate {
    my $p = shift;

    my $model = $p->{'model'};
    my $plate_name_input = $p->{'plate_name'};
    my $repeat_mask = $p->{'repeat_mask_list'};
    my $species_input = $p->{'species'};
    my $assembly_input = $p->{'assembly'};
    my $plate_well_input = $p->{'plate_well'};
    my $crispr_pair_id = $p->{'crispr_pair_id'};
    my $plate_crispr_right = $p->{'plate_crispr_right'};
    my $plate_crispr_left = $p->{'plate_crispr_left'};

    $logger->info("Starting primer generation for plate $plate_name_input");
    my $rpt_string = join( ',', @$repeat_mask);
    $logger->info("Using sequence repeat mask of: $rpt_string");


    # Probably get the design ids from a file or the command line or directly from the database.
    my $species = $species_input;

    my $plate_rs = $model->schema->resultset( 'Plate' )->search(
        {
            'name' => $plate_name_input
        },
    );

     if ( !$plate_rs->count ) {
         $logger->fatal("No wells on plate or non-existent plate: $plate_name_input");
         exit;
     }


    my $plate = $plate_rs->first;
    my $plate_name = $plate->name;
    $logger->info( 'Plate selected: ' . $plate_name );

    my @wells = $plate->wells->all;

    my $well_count = @wells;
    $logger->info( 'Processing crispr primers for ' . $well_count . ' wells:');
    my @well_id_list;

    # Select just the well we want
    if ( $plate_well_input ) {
        $logger->info("Selecting well: [$plate_well_input]");
        @wells = grep { $_->name =~ /$plate_well_input/  } @wells;
    }

    foreach my $well ( @wells ) {
        push @well_id_list, $well->id;
    }

    my $design_data_cache = $model->create_design_data_cache(
            \@well_id_list,
        );

    $logger->debug( 'Created design well data cache' );



    my $dc = run_primers({
            'wells' => \@wells,
            'design_data_cache' => $design_data_cache,
            'model' => $model,
            'species' => $species,
            'plate_name' => $plate_name,
            'repeat_mask' => $repeat_mask,
            'assembly_id' => $assembly_input,
            'plate_well' => $plate_well_input,
            'crispr_pair_id' => $crispr_pair_id,
            'plate_crispr_right' => $plate_crispr_right,
            'plate_crispr_left' => $plate_crispr_left,
        });

    $logger->info( 'Done' );
    return;
}

#

sub create_d_plate_hash{
    my $model = shift;
    my $d_plate_list = shift;

    my $dis_designs;
    my @well_id_list;

    foreach my $d_plate_name ( @{$d_plate_list} ) {
        my $plate_rs = $model->schema->resultset( 'Plate' )->search(
            {
                'name' => $d_plate_name
            },
        );

        if ( !$plate_rs->count ) {
            $logger->fatal("No wells on disambiguation plate or non-existent plate: $d_plate_name");
            die;
        }


        my $plate = $plate_rs->first;
        my $plate_name = $plate->name;
        $logger->info( 'Disambiguation plate selected: ' . $plate_name );

        my @wells = $plate->wells->all;


        foreach my $well ( @wells ) {
            push @well_id_list, $well->id;
        }
    }

    $dis_designs = $model->create_design_data_cache(
                \@well_id_list,
            );

    return $dis_designs;
}

sub verify_and_update_design_data_cache{
    my $well_id = shift;
    my $crispr_id = shift;
    my $dd_cache = shift;
    my $dis_designs = shift;
    my $schema = shift;
    my $assembly = shift;

    # In a special edge case, the design data cache will not be populated
    # because designs were not added to the plate. We need to use the crispr locus data for each well
    # to backtrack and locate the design so that the design data cache can be populated.
    # Then we can process as normal.
    my $changed = 0;

    if ( ! exists $dd_cache->{$well_id}->{'design_id'} ){
        $logger->debug( 'No design id value found for well_id: ' . $well_id);
        # get the design for the crispr
        if ( my ($design_id, $gene_id ) = get_design_for_single_crispr( $crispr_id, $dis_designs, $schema, $assembly ) ){
           $dd_cache->{$well_id}->{'design_id'} = $design_id;
           $dd_cache->{$well_id}->{'gene_id'} = $gene_id;
           $changed = 1;
        }
        else {
            $logger->info( 'No design information found for well_id: ' . $well_id . 'crispr_id: ' . $crispr_id);
        }
    }

    return $dd_cache, $changed;
}

sub get_design_for_single_crispr {
    my $crispr_id = shift;
    my $dis_designs = shift;
    my $schema = shift;
    my $assembly = shift;

    my %crispr_data;
    my $crispr = single_crispr_oligo_result( $schema, $crispr_id );

    my $locus_count = $crispr->loci->count;
    if ($locus_count != 1 ) {
        INFO ('Found multiple loci for ' . $crispr_id);
        INFO ('Using first found locus' );
    }
    my $locus = $crispr->loci->first;
    $crispr_data{'chr_start'} = $locus->chr_start;
    $crispr_data{'chr_end'} = $locus->chr_end;
    # Do not use the crispr strand information! Use it from the design for the gene
    $crispr_data{'chr_id'} = $locus->chr_id;
    $crispr_data{'chr_name'} = $locus->chr->name;
    $crispr_data{'seq'} = $crispr->seq;


    my $design_rs = $schema->resultset('GenericDesignBrowser')->search( {},
        {
            bind => [
                $crispr_data{'chr_start'} - 2000,
                $crispr_data{'chr_end'} + 2000,
                $crispr_data{'chr_id'},
                $assembly,
            ],
        }
    );
    my $design_data;
    my @d_dis_well_keys = keys %{$dis_designs};
    my @d_dis_design_ids;
    foreach my $well_key ( @d_dis_well_keys ) {
        push @d_dis_design_ids, $dis_designs->{$well_key}->{'design_id'};
    }

    my @design_results;
    if ($design_rs->count > 1) {
        $logger->warn('Multiple designs found for crispr_id: ' . $crispr_id);
        my $counter = 0;
        while ( my $row = $design_rs->next ) {
            $counter++;
            $logger->warn($counter . ': Design id: ' . $row->design_id);
            # check which is the correct design in the disambiguation design hash
            my @design_matches = grep {
                $_ == $row->design_id
            } @d_dis_design_ids;
            if ( scalar @design_matches < 1 ) {
                $logger->warn('No disambiguation designs match');
            }
            else {
                push @design_results, $row;
            }
        }
        # If there was only one design row use it, otherwise disambiguation failed.
        if (scalar @design_results > 1 ){
            $logger->error('Disambiguation failed..');
            foreach my $design_r ( @design_results ) {
                $logger->warn('Potential Design: ' . $design_r->design_id);
            }
            # FIXME: temp fudge to keep things going.
            $design_data = $design_results[0];
        }
        else {
            $design_data = $design_results[0];
        }
    }
    else {
        # There was only one design
        $design_data = $design_rs->first;
    }
    my $design_id = $design_data->design_id;
    my $gene_id = $design_data->gene_id;

    return $design_id, $gene_id
}

sub single_crispr_oligo_result {
    my $schema = shift;
    my $crispr_id = shift;

    my $crispr_rs = $schema->resultset('Crispr')->find(
        {
            'id' => $crispr_id,
        },
    );

    return $crispr_rs;
}

sub run_single_crispr_primers {
    my $params = shift;
    my $wells = $params->{'wells'};
    my $design_data_cache = $params->{'design_data_cache'};
    my $model = $params->{'model'};
    my $species = $params->{'species'};
    my $assembly = $params->{'assembly'};
    my $plate_name = $params->{'plate_name'};
    my $repeat_mask = $params->{'repeat_mask'};
    my $dis_designs = $params->{'dis_designs'};

    my $formatted_well = $params->{'plate_well'} ? ('_' . $params->{'plate_well'}) : '';

    my $lines;
    $logger->debug( 'Generating gene symbol cache' );

    $design_data_cache = generate_gene_symbols_cache({
            'wells' => $wells,
            'design_data_cache' => $design_data_cache,
            'species' => $species,
            'model' => $model,
        });
     my ($out_rows, $crispr_clip );
     if ( uc($crispr_primers_required) eq 'YES' ) {
        $logger->debug( 'Preparing crispr primers');
        ($out_rows, $crispr_clip ) = prepare_single_crispr_primers({
                'model' => $model,
                'wells' => $wells,
                'design_data_cache' => $design_data_cache,
                'species' => $species,
                'assembly' => $assembly,
                'repeat_mask' => $repeat_mask,
                'dis_designs' => $dis_designs,
            });
        $logger->debug( 'Generating single crispr primer output file' );
        $lines = generate_single_crispr_output( $out_rows );
        create_output_file( $plate_name . $formatted_well . '_single_crispr_primers.csv', $lines );
    }
    else {
        $logger->info( 'Crispr primers not required' );
    }

    if ( uc($pcr_primers_required) eq 'YES' ) {
        $logger->info( 'Generating PCR primers for crispr region' );
        ($out_rows, $crispr_clip) = prepare_pcr_primers({
                'model' => $model,
                'crispr_clip' => $crispr_clip,
                'wells' => $wells,
                'design_data_cache' => $design_data_cache,
                'species' => $species,
                'repeat_mask' => $repeat_mask,
                'assembly_id' => $assembly,
            });
        $logger->debug( 'Generating pcr primer output file' );
        $lines = generate_pcr_output( $out_rows );
        create_output_file( $plate_name . $formatted_well . '_pcr_primers.csv', $lines );
    }
    else {
        $logger->info( 'PCR primers not required' );
    }

    if ( uc($genotyping_primers_required) eq 'YES' ) {
        $logger->info( 'Generating genotyping primers' );
        my $design_oligos;
        ($out_rows, $crispr_clip) = prepare_genotyping_primers({
                'model' => $model,
                'crispr_clip' => $crispr_clip,
                'wells' => $wells,
                'design_data_cache' => $design_data_cache,
                'species' => $species,
                'repeat_mask' => $repeat_mask,
                'assembly_id' => $assembly,
            });

        $lines = generate_genotyping_output( $out_rows );
        $logger->info( 'Generating genotyping primer output file' );
        create_output_file( $plate_name . $formatted_well . '_genotyping_primers.csv' ,$lines );
    }
    else {
        $logger->info( 'Genotyping primers not required' );
    }

    return $crispr_clip;
}

sub run_primers {
    my $params = shift;

    my $wells = $params->{'wells'};
    my $design_data_cache = $params->{'design_data_cache'};
    my $model = $params->{'model'};
    my $species = $params->{'species'};
    my $plate_name = $params->{'plate_name'};
    my $plate_crispr_left = $params->{'plate_crispr_left'};
    my $plate_crispr_right = $params->{'plate_crispr_right'};
    my $repeat_mask= $params->{'repeat_mask'};
    my $crispr_pair_id_inp = $params->{'crispr_pair_id'};
    my $assembly_id = $params->{'assembly_id'};

    my $formatted_well = $params->{'plate_well'} ? ('_' . $params->{'plate_well'}) : '';

    my $lines;
    $logger->debug( 'Generating gene symbol cache' );

    $design_data_cache = generate_gene_symbols_cache({
            'wells' => $wells,
            'design_data_cache' => $design_data_cache,
            'species' => $species,
            'model' => $model,
        });
    my ($out_rows, $crispr_clip ); 
    if ( uc($crispr_primers_required) eq 'YES') {
        $logger->debug( 'Preparing crispr primers');
        ($out_rows, $crispr_clip ) = prepare_crispr_primers({
                'model' => $model,
                'wells' => $wells,
                'design_data_cache' => $design_data_cache,
                'species' => $species,
                'repeat_mask' => $repeat_mask,
                'crispr_pair_id' => $crispr_pair_id_inp,
                'assembly_id' => $assembly_id,
                'plate_crispr_left' => $plate_crispr_left,
                'plate_crispr_right' => $plate_crispr_right,
            });
        $logger->debug( 'Generating crispr primer output file' );
        $lines = generate_crispr_output( $out_rows );
        create_output_file( $plate_name . $formatted_well . '_crispr_primers.csv', $lines );
    }
    else {
        $logger->info( 'Crispr primers not required' );
    }

    if ( uc($pcr_primers_required) eq 'YES') {
        $logger->info( 'Generating PCR primers for crispr region' );
        ($out_rows, $crispr_clip) = prepare_pcr_primers({
                'model' => $model,
                'crispr_clip' => $crispr_clip,
                'wells' => $wells,
                'design_data_cache' => $design_data_cache,
                'species' => $species,
                'repeat_mask' => $repeat_mask,
                'assembly_id' => $assembly_id,
            });
        $logger->debug( 'Generating pcr primer output file' );
        $lines = generate_pcr_output( $out_rows );
        create_output_file( $plate_name . $formatted_well . '_pcr_primers.csv', $lines );
    }
    else {
        $logger->info( 'PCR primers not required' );
    }

    if ( uc($genotyping_primers_required) eq 'YES') {
        $logger->info( 'Generating genotyping primers' );
        my $design_oligos;
        ($out_rows, $crispr_clip) = prepare_genotyping_primers({
                'model' => $model,
                'crispr_clip' => $crispr_clip,
                'wells' => $wells,
                'design_data_cache' => $design_data_cache,
                'species' => $species,
                'repeat_mask' => $repeat_mask,
                'assembly_id' => $assembly_id,
            });

        $lines = generate_genotyping_output( $out_rows );
        $logger->info( 'Generating genotyping primer output file' );
        create_output_file( $plate_name . $formatted_well . '_genotyping_primers.csv' ,$lines );
    }
    else {
        $logger->info( 'Genotyping primers not required' );
    }

    return $crispr_clip;
}


sub prepare_pcr_primers {
    my $params = shift;

    my $model = $params->{'model'};
    my $crispr_clip = $params->{'crispr_clip'};
    my $wells = $params->{'wells'};
    my $design_data_cache = $params->{'design_data_cache'};
    my $species = $params->{'species'};
    my $repeat_mask = $params->{'repeat_mask'};
    my $assembly_id = $params->{'assembly_id'};

    foreach my $well ( @{$wells} ) {

        my $well_id = $well->id;
        my $well_name = $well->name;

        my $gene_name = $design_data_cache->{$well_id}->{'gene_symbol'};

        $logger->info( 'crispr_pcr: ' . $well_name . "\t(" . $gene_name . ')' );
        next if $crispr_clip->{$well_name}->{'crispr_primers'}->{'error_flag'} ne 'pass';

        my ($crispr_pcr_primers, $crispr_pcr_mapped, $chr_seq_start)
            = LIMS2::Model::Util::OligoSelection::pick_crispr_PCR_primers( {
                schema => $model->schema,
                well_id => $well,
                crispr_primers => $crispr_clip->{$well_name},
                species => $species,
                repeat_mask => $repeat_mask,
            });
        $crispr_clip->{$well_name}->{'crispr_pcr_primers'} = $crispr_pcr_mapped;
        $crispr_clip->{$well_name}{'chr_seq_start'} = $chr_seq_start;

    }
    my @out_rows;
    my $csv_row;
    my $primer_type = 'crispr_pcr_primers';
    foreach my $well_name ( keys %{$crispr_clip} ) {
        my @out_vals = (
            $well_name,
            $crispr_clip->{$well_name}->{'gene_name'} // 'not set',
            $crispr_clip->{$well_name}->{'design_id'} // 'not set',
            $crispr_clip->{$well_name}->{'strand'} // 'not set',
        );
        if ( $crispr_clip->{$well_name}->{'crispr_primers'}->{'error_flag'} ne 'pass' ){

            $csv_row = join( ',' , @out_vals);
            push @out_rows, $csv_row;
            
            next;
        }
        # Take the two highest ranking primers
        my ($rank_a, $rank_b) = get_best_two_primer_ranks( $crispr_clip->{$well_name}->{$primer_type}->{'left'} );
        # only need left or right as both will be the same for rank purposes
        if ($crispr_clip->{$well_name}->{'strand'} eq 'plus') {
            # first set
            push (@out_vals,
                data_to_push($crispr_clip, $primer_type, $well_name, 'left', $rank_a // '99')
            );
            persist_pcr_primers($model, 'PF1', $crispr_clip, $primer_type, $well_name, 'left', $rank_a // '99', $assembly_id);
            push (@out_vals, (
                data_to_push($crispr_clip, $primer_type, $well_name, 'right', $rank_a // '99')
            ));
            persist_pcr_primers($model, 'PR1', $crispr_clip, $primer_type, $well_name, 'right', $rank_a // '99', $assembly_id);
            # second set
            push (@out_vals,
                data_to_push($crispr_clip, $primer_type, $well_name, 'left', $rank_b // '99')
            );
            persist_pcr_primers($model, 'PF2', $crispr_clip, $primer_type, $well_name, 'left', $rank_b // '99', $assembly_id);
            push (@out_vals, (
                data_to_push($crispr_clip, $primer_type, $well_name, 'right', $rank_b // '99')
            ));
            persist_pcr_primers($model, 'PR2', $crispr_clip, $primer_type, $well_name, 'right', $rank_b // '99', $assembly_id);
        }
        elsif ($crispr_clip->{$well_name}->{'strand'} eq 'minus') {
            push (@out_vals,
                data_to_push($crispr_clip, $primer_type, $well_name, 'right', $rank_a // '99')
            );
            persist_pcr_primers($model, 'PF1', $crispr_clip, $primer_type, $well_name, 'right', $rank_a // '99', $assembly_id);
            push (@out_vals, (
                data_to_push($crispr_clip, $primer_type, $well_name, 'left', $rank_a // '99')
            ));
            persist_pcr_primers($model, 'PR1', $crispr_clip, $primer_type, $well_name, 'left', $rank_a // '99', $assembly_id);
            # PF2/PR2
            push (@out_vals,
                data_to_push($crispr_clip, $primer_type, $well_name, 'right', $rank_b // '99')
            );
            persist_pcr_primers($model, 'PF2', $crispr_clip, $primer_type, $well_name, 'right', $rank_b // '99', $assembly_id);
            push (@out_vals, (
                data_to_push($crispr_clip, $primer_type, $well_name, 'left', $rank_b // '99')
            ));
            persist_pcr_primers($model, 'PR2', $crispr_clip, $primer_type, $well_name, 'left', $rank_b // '99', $assembly_id);

        }
        else {
            $logger->error( 'Chromosome strand not defined.' );
        }
        # Anything after here for context only

        $csv_row = join( ',' , @out_vals);
        push @out_rows, $csv_row;
    }
    return \@out_rows, $crispr_clip;
}

sub prepare_genotyping_primers {
    my $params = shift;

    my $model = $params->{'model'};
    my $primer_clip = $params->{'crispr_clip'};
    my $wells = $params->{'wells'};
    my $design_data_cache = $params->{'design_data_cache'};
    my $species = $params->{'species'};
    my $repeat_mask = $params->{'repeat_mask'};
    my $assembly_id = $params->{'assembly_id'};


#    my %primer_clip;

    my $design_row;
    foreach my $well ( @{$wells} ) {
        my $well_id = $well->id;
        my $design_id = $design_data_cache->{$well_id}->{'design_id'};
        my $gene_id = $design_data_cache->{$well_id}->{'gene_id'};
        my $well_name = $well->name;
        my $gene_name;

        $gene_name = $design_data_cache->{$well_id}->{'gene_symbol'};

        $logger->info( $design_id . "\t(" . $gene_name . ')' );
        next if $primer_clip->{$well_name}->{'crispr_primers'}->{'error_flag'} ne 'pass';

        my ($genotyping_primers, $genotyping_mapped, $chr_strand, $design_oligos, $chr_seq_start)
            = LIMS2::Model::Util::OligoSelection::pick_genotyping_primers( {
                schema => $model->schema,
                design_id => $design_id,
                well_id => $well,
                species => $species,
                repeat_mask => $repeat_mask,
            });
        # TODO: check whether already exists first? Should have already been filled in?
        $primer_clip->{$well_name}{'gene_name'} = $gene_name;
        $primer_clip->{$well_name}{'design_id'} = $design_id;
        $primer_clip->{$well_name}{'strand'} = $chr_strand;

        $primer_clip->{$well_name}{'genotyping_primers'} = $genotyping_mapped;
        $primer_clip->{$well_name}{'design_oligos'} = $design_oligos;
        $primer_clip->{$well_name}{'chr_seq_start'} = $chr_seq_start; # This changes for each call
    }
    my @out_rows;
    my $csv_row;
    my $primer_type = 'genotyping_primers';
    foreach my $well_name ( keys %$primer_clip ) {
        my @out_vals = (
            $well_name,
            $primer_clip->{$well_name}->{'gene_name'},
            $primer_clip->{$well_name}->{'design_id'},
            $primer_clip->{$well_name}->{'strand'},
        );
        if ( $primer_clip->{$well_name}->{'crispr_primers'}->{'error_flag'} ne 'pass' ){

            $csv_row = join( ',' , @out_vals);
            push @out_rows, $csv_row;
            
            next;
        }
        # Take the two highest ranking primers
        my ($rank_a, $rank_b) = get_best_two_primer_ranks( $primer_clip->{$well_name}->{$primer_type}->{'left'} );
        # only need left or right as both will be the same for rank purposes
        # GF1/GR1
        if ($primer_clip->{$well_name}->{'strand'} eq 'plus') {
            push (@out_vals,
                data_to_push($primer_clip, $primer_type, $well_name, 'left', $rank_a // '99')
            );
            persist_gen_primers($model, 'GF1', $primer_clip, $primer_type, $well_name, 'left', $rank_a // '99', $assembly_id);
            push (@out_vals, (
                data_to_push($primer_clip, $primer_type, $well_name, 'right', $rank_a // '99')
            ));
            persist_gen_primers($model, 'GR1', $primer_clip, $primer_type, $well_name, 'right', $rank_a // '99', $assembly_id);
            # GF2/GR2
            push (@out_vals,
                data_to_push($primer_clip, $primer_type, $well_name, 'left', $rank_b // '99')
            );
            persist_gen_primers($model, 'GF2', $primer_clip, $primer_type, $well_name, 'left', $rank_b // '99', $assembly_id);
            push (@out_vals, (
                data_to_push($primer_clip, $primer_type, $well_name, 'right', $rank_b // '99')
            ));
            persist_gen_primers($model, 'GR2', $primer_clip, $primer_type, $well_name, 'right', $rank_b // '99', $assembly_id);
        }
        elsif ($primer_clip->{$well_name}->{'strand'} eq 'minus') {
            push (@out_vals,
                data_to_push($primer_clip, $primer_type, $well_name, 'right', $rank_a // '99')
            );
            persist_gen_primers($model, 'GF1', $primer_clip, $primer_type, $well_name, 'right', $rank_a // '99', $assembly_id);
            push (@out_vals, (
                data_to_push($primer_clip, $primer_type, $well_name, 'left', $rank_a // '99')
            ));
            persist_gen_primers($model, 'GR1', $primer_clip, $primer_type, $well_name, 'left', $rank_a // '99', $assembly_id);
            # GF2/GR2
            push (@out_vals,
                data_to_push($primer_clip, $primer_type, $well_name, 'right', $rank_b // '99')
            );
            persist_gen_primers($model, 'GF2', $primer_clip, $primer_type, $well_name, 'right', $rank_b // '99', $assembly_id);
            push (@out_vals, (
                data_to_push($primer_clip, $primer_type, $well_name, 'left', $rank_b // '99')
            ));
            persist_gen_primers($model, 'GR2', $primer_clip, $primer_type, $well_name, 'left', $rank_b // '99', $assembly_id);
        }
        else {
            $logger->error( 'Chromosome strand not defined.' );
        }
        my @d_oligo_keys = sort keys %{$primer_clip->{$well_name}{'design_oligos'}};
        push (@out_vals,
            $d_oligo_keys[0],
            $primer_clip->{$well_name}->{'design_oligos'}->{$d_oligo_keys[0]}->{'seq'},
            $d_oligo_keys[1],
            $primer_clip->{$well_name}->{'design_oligos'}->{$d_oligo_keys[1]}->{'seq'},
        );

        $csv_row = join( ',' , @out_vals);
        push @out_rows, $csv_row;
    }

    return \@out_rows, $primer_clip;
}

sub data_to_push {
    my $pc = shift;
    my $primer_type = shift;
    my $well_name = shift;
    my $lr = shift;
    my $rank = shift;

    my $primer_name = $lr . '_' . $rank;
    return (
        $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'seq'} // "'-",
        $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'location'}->{'_strand'} // "'-",
        $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'length'} // "'-",
        $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'gc_content'} // "'-",
        $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'melting_temp'} // "'-",
    );
}

=head persist_primers

Write the calculated data to the store

=cut

sub persist_gen_primers {
    return persist_primers( @_, 'genotyping' );
}

sub persist_seq_primers {
    return persist_primers( @_, 'sequencing');
}

sub persist_pcr_primers {
    return persist_primers( @_, 'pcr');
}

sub persist_primers {
    #Global param check
    return if ( $persist_param ne 'DB' && $persist_param ne 'BOTH' );
    my $model = shift;
    my $primer_label = shift;

    my $pc = shift;
    my $primer_type = shift;
    my $well_name = shift;
    my $lr = shift;
    my $rank = shift;
    my $assembly_id = shift;
    my $primer_class = shift;

    return if ( $rank eq '99'); # no suitable primers generated
    $logger->info("Persisting $primer_class primers for well $well_name ...");
    my $primer_name = $lr . '_' . $rank;
    my $seq = $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'seq'};
    my $gc = $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'gc_content'};
    my $tm = $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'melting_temp'};
    my $chr_strand = $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'location'}->{'_strand'};
    my $chr_start = $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'location'}->{'_start'}
        + $pc->{$well_name}->{'chr_seq_start'} - 1;
    my $chr_end = $pc->{$well_name}->{$primer_type}->{$lr}->{$primer_name}->{'location'}->{'_end'}
        + $pc->{$well_name}->{'chr_seq_start'} - 1;
    my $chr_id = $pc->{$well_name}->{'crispr_seq'}->{'left_crispr'}->{'chr_id'}; # not the translated name

    my $crispr_primer_result;
    if ($primer_class eq 'genotyping' ) {
        if ( $seq ) {
            # There is no unique constraint on the GenotypingPrimer resultset.
            # Therefore, check whether this combination already exists before updating
            # the table
            # This must be done in a transaction to prevent a race condition
            my $coderef = sub {
                my $crispr_primer_result_set =
                    $model->schema->resultset( 'GenotypingPrimer' )->search({
                            'design_id' => $pc->{$well_name}->{design_id},
                            'genotyping_primer_type_id' => $primer_label,
                        });
                if ($crispr_primer_result_set->count == 0 ) {
                    $crispr_primer_result = $model->schema->resultset('GenotypingPrimer')->create({
                    'design_id'      => $pc->{$well_name}->{'design_id'},
                    'genotyping_primer_type_id'
                                     => $primer_label,
                    'seq'            => $seq,
                    'tm'             => $tm ,
                    'gc_content'     => $gc,
                        'genotyping_primer_loci' => [{
                             "assembly_id" => $assembly_id,
                             "chr_id" => $chr_id,
                             "chr_start" => $chr_start,
                             "chr_end" => $chr_end,
                             "chr_strand" => $chr_strand,
                         }],
                    });
                }
                else {
                    $logger->warn('Genotyping design/primer already exists: '
                        . $pc->{$well_name}->{design_id}
                        . '/'
                        . $primer_label
                    );
                }
                return $crispr_primer_result;
            };
            #
            my $rs;
            try {
                $rs = $model->schema->txn_do( $coderef );
            } catch {
                my $error = shift;
                # Transaction failed
                die 'Something went wrong with that transaction: ' . $error;
            };
        }

    }
    else { # it must be sequencing or pcr
        if ( $seq ) {
            my $create_params = {
                'primer_name'    => $primer_label,
                'primer_seq'     => $seq,
                'tm'             => $tm ,
                'gc_content'     => $gc,
                    'crispr_primer_loci' => [{
                         "assembly_id" => $assembly_id,
                         "chr_id" => $chr_id,
                         "chr_start" => $chr_start,
                         "chr_end" => $chr_end,
                         "chr_strand" => $chr_strand,
                     }],
            };
            my $search_params = ();
            if ( $pc->{$well_name}->{'pair_id'} ) {
                $create_params->{'crispr_pair_id'} = $pc->{$well_name}->{'pair_id'};
                $search_params->{'crispr_pair_id'} = $pc->{$well_name}->{'pair_id'};
            }
            else {
                $create_params->{'crispr_id'} = $pc->{$well_name}->{'crispr_id'};
                $search_params->{'crispr_id'} = $pc->{$well_name}->{'crispr_id'};
            }
            # Can't use update_or_create, even though this is a 1-1 relationship, dbix thinks it is multi.
            # Therefore need to check whether exists and update it myself if it does, otherwise create a new
            # database tuple.
            $search_params->{'primer_name'} = $primer_label;
            my $coderef = sub {
                my $crispr_check_r = $model->schema->resultset('CrisprPrimer')->find( $search_params );
                if ( $crispr_check_r ) {
                    $logger->info( 'Deleting entry for '
                        . ($search_params->{'crispr_pair_id'} // $search_params->{'crispr_id'} )
                        . ' label: '
                        . $search_params->{'primer_name'}
                        . ' from the database');
                    if ( ! $crispr_check_r->delete ) {
                        $logger->info( 'Unable to delete record(s) from the database' );
                        die;
                    }
                }
                $crispr_primer_result = $model->schema->resultset('CrisprPrimer')->create( $create_params );
                if ( ! $crispr_primer_result ) {
                    $logger->info('Unable to create crispr primer records for '
                        . ($search_params->{'crispr_pair_id'} // $search_params->{'crispr_id'})
                        . ' label: '
                        . $search_params->{'primer_name'}
                    );
                }
                else {
                    $logger->info('Created '
                        . ($search_params->{'crispr_pair_id'} // $search_params->{'crispr_id'})
                        . ' label: '
                        . $search_params->{'primer_name'}
                    );
                }
                return $crispr_primer_result;
            };
            my $rs;
            try {
                $rs = $model->schema->txn_do( $coderef );
            } catch {
                my $error = shift;
                # Transaction failed
                die 'Something went wrong with that transaction: ' . $error;
            };
        }
    } # End if genotyping
#  print $crispr_primer->in_storage();  ## 1 (TRUE)
    return $crispr_primer_result;
}

sub get_best_two_primer_ranks {
    my $pc_rank = shift;
    my @rank_keys = sort keys %$pc_rank;

    s/left_// for @rank_keys;
    return ($rank_keys[0], $rank_keys[1] );
}

=head prepare_single_crispr_primers

=cut

sub prepare_single_crispr_primers {
    my $params = shift;

    my $model = $params->{'model'};
    my $wells = $params->{'wells'};
    my $design_data_cache = $params->{'design_data_cache'};
    my $species = $params->{'species'};
    my $assembly = $params->{'assembly'};
    my $repeat_mask = $params->{'repeat_mask'};
    my $dis_designs = $params->{'dis_designs'};

    my %primer_clip;

    my $design_row;
    foreach my $well ( @{$wells} ) {
        my $well_id = $well->id;
        my $well_name = $well->name;
        my $gene_name;

        my ($crispr_left, $crispr_right) = $well->left_and_right_crispr_wells;
        my $crispr_id = $crispr_left->crispr->id;
        ($design_data_cache, my $changed) = verify_and_update_design_data_cache( $well_id, $crispr_id, $design_data_cache, $dis_designs, $model->schema, $assembly );
        if ($changed) {
            # update gene symbols
            my @wlist = ( $well );
            $design_data_cache = generate_gene_symbols_cache({
                'wells' => \@wlist,
                'design_data_cache' => $design_data_cache,
                'species' => $species,
                'model' => $model,
                });
        }
        my $design_id = $design_data_cache->{$well_id}->{'design_id'};
        my $gene_id = $design_data_cache->{$well_id}->{'gene_id'};
        $gene_name = $design_data_cache->{$well_id}->{'gene_symbol'};
        $logger->info( "$design_id\t$gene_name\tcrispr_id:\t$crispr_id" );

        my ($crispr_results, $crispr_primers, $chr_strand) = LIMS2::Model::Util::OligoSelection::pick_single_crispr_primers( {
                schema => $model->schema,
                design_id => $design_id,
                crispr_id => $crispr_id,
                species => $species,
                repeat_mask => $repeat_mask,
            });
        $primer_clip{$well_name}{'crispr_id'} = $crispr_id;
        $primer_clip{$well_name}{'gene_name'} = $gene_name;
        $primer_clip{$well_name}{'design_id'} = $design_id;
        $primer_clip{$well_name}{'strand'} = $chr_strand;

        $primer_clip{$well_name}{'crispr_seq'} = $crispr_results;
        $primer_clip{$well_name}{'crispr_primers'} = $crispr_primers;
    }

    my @out_rows;
    my $rank = 0; # always show the rank 0 primer
    my $primer_type = 'crispr_primers';
    foreach my $well_name ( keys %primer_clip ) {
        my @out_vals = (
            $well_name,
            $primer_clip{$well_name}{'gene_name'},
            $primer_clip{$well_name}{'design_id'},
            $primer_clip{$well_name}{'strand'},
            $primer_clip{$well_name}{'crispr_id'},
        );
        push (@out_vals, (
            data_to_push(\%primer_clip, $primer_type, $well_name, 'left', $rank // '99')
        ));
        persist_seq_primers($model, 'SF1', \%primer_clip, $primer_type, $well_name, 'left', $rank // '99', $assembly);
        push (@out_vals, (
            data_to_push(\%primer_clip, $primer_type, $well_name, 'right', $rank // '99')
        ));
        persist_seq_primers($model, 'SR1', \%primer_clip, $primer_type, $well_name, 'right', $rank // '99', $assembly);
        push (@out_vals,
            $primer_clip{$well_name}->{'crispr_seq'}->{'left_crispr'}->{'id'},
            $primer_clip{$well_name}->{'crispr_seq'}->{'left_crispr'}->{'seq'},
        );
        my $csv_row = join( ',' , @out_vals);
        push @out_rows, $csv_row;
    }

    return (\@out_rows, \%primer_clip);
}

sub prepare_crispr_primers {
    my $params = shift;

    my $model = $params->{'model'};
    my $wells = $params->{'wells'};
    my $design_data_cache = $params->{'design_data_cache'};
    my $species = $params->{'species'};
    my $repeat_mask = $params->{'repeat_mask'};
    my $crispr_pair_id_inp = $params->{'crispr_pair_id'} // '';
    my $assembly_id = $params->{'assembly_id'};
    my $plate_crispr_left = $params->{'plate_crispr_left'};
    my $plate_crispr_right = $params->{'plate_crispr_right'};
    my %primer_clip;

    my $crispr_pair_id_cache = generate_crispr_pair_id_cache($model, $plate_crispr_left, $plate_crispr_right);

    my $design_row;
    foreach my $well ( @{$wells} ) {
        my $well_id = $well->id;
        my $design_id = $design_data_cache->{$well_id}->{'design_id'};
        my $gene_id = $design_data_cache->{$well_id}->{'gene_id'};
        my $well_name = $well->name;
        my $gene_name;

        $gene_name = $design_data_cache->{$well_id}->{'gene_symbol'};

        my $crispr_pair_id;
        if ( $crispr_pair_id_inp eq '' ) {
            # Get the pair information by resolving the left and right crispr plate information
            if ( $crispr_pair_id_cache->{$well_name} ) {
                $crispr_pair_id = $crispr_pair_id_cache->{$well_name}->{'crispr_pair_id'};
            }
            # Get the crispr_pair and thence its id from the well
            elsif ( my $crispr_pair = $well->crispr_pair ) {
                $crispr_pair_id = $crispr_pair->id;
            }
        }
        else {
            $crispr_pair_id = $crispr_pair_id_inp;
        }

        $logger->info( "$design_id\t$gene_name\tcrispr_pair_id:\t$crispr_pair_id" );
        # Check whether there are already primers available for this crispr_pair_id
        if ( check_primers_for_crispr_pair_id( $model, $crispr_pair_id )->count > 0 ) {
            $logger->warn( '++++++++ primers already exist for crispr_pair_id: ' . $crispr_pair_id );
            # TODO: Add code to alter persistence behaviour
        }

        my ($crispr_results, $crispr_primers, $chr_strand, $chr_seq_start) = LIMS2::Model::Util::OligoSelection::pick_crispr_primers( {
                'schema' => $model->schema,
                'design_id' => $design_id,
                'crispr_pair_id' => $crispr_pair_id,
                'species' => $species,
                'repeat_mask' => $repeat_mask,
            });
        $primer_clip{$well_name}{'pair_id'} = $crispr_pair_id;
        $primer_clip{$well_name}{'gene_name'} = $gene_name;
        $primer_clip{$well_name}{'design_id'} = $design_id;
        $primer_clip{$well_name}{'strand'} = $chr_strand;

        $primer_clip{$well_name}{'crispr_seq'} = $crispr_results;
        $primer_clip{$well_name}{'crispr_primers'} = $crispr_primers;
        $primer_clip{$well_name}{'chr_seq_start'} = $chr_seq_start;
    }

    my @out_rows;
    my $csv_row;
    my $rank = 0; # always show the rank 0 primer
    my $primer_type = 'crispr_primers';
    foreach my $well_name ( keys %primer_clip ) {
        my @out_vals = (
            $well_name,
            $primer_clip{$well_name}{'gene_name'},
            $primer_clip{$well_name}{'design_id'},
            $primer_clip{$well_name}{'strand'},
            $primer_clip{$well_name}{'pair_id'},
        );
        if ( $primer_clip{$well_name}->{'crispr_primers'}->{'error_flag'} ne 'pass' ){
            push @out_vals
                , 'primer3_explain_left'
                , $primer_clip{$well_name}->{'crispr_primers'}->{'primer3_explain_left'}
                , 'primer3_explain_right'
                , $primer_clip{$well_name}->{'crispr_primers'}->{'primer3_explain_right'};

            $csv_row = join( ',' , @out_vals);
            push @out_rows, $csv_row;
            
            next;
        }

        push (@out_vals, (
            data_to_push(\%primer_clip, $primer_type, $well_name, 'left', $rank // '99')
        ));
        persist_seq_primers($model, 'SF1', \%primer_clip, $primer_type, $well_name, 'left', $rank // '99', $assembly_id);
        push (@out_vals, (
            data_to_push(\%primer_clip, $primer_type, $well_name, 'right', $rank // '99')
        ));
        persist_seq_primers($model, 'SR1', \%primer_clip, $primer_type, $well_name, 'right', $rank // '99', $assembly_id);
        push (@out_vals,
            $primer_clip{$well_name}->{'crispr_seq'}->{'left_crispr'}->{'id'},
            $primer_clip{$well_name}->{'crispr_seq'}->{'left_crispr'}->{'seq'},
        );
        push (@out_vals,
            $primer_clip{$well_name}->{'crispr_seq'}->{'right_crispr'}->{'id'},
            $primer_clip{$well_name}->{'crispr_seq'}->{'right_crispr'}->{'seq'},
        );
        my $csv_row = join( ',' , @out_vals);
        push @out_rows, $csv_row;
    }

    return (\@out_rows, \%primer_clip);
}


=head check_primers_for_crispr_pair_id( $model, $crispr_pair_id )
Given:
    model,
    crispr_pair_id

Returns:
    resultset

Look up db entry for primers for this crispr_pair_id

=cut

sub check_primers_for_crispr_pair_id {
    my $model = shift;
    my $crispr_pair_id = shift;

    my $crispr_primers_rs = $model->schema->resultset('CrisprPrimer')->search({
        'crispr_pair_id' => $crispr_pair_id,
    });

    return $crispr_primers_rs;
}

=head generate_fsa_file
Generates data to output fasta format files for use with BlastN

=cut

sub generate_fsa_files {
    my $clip = shift;

    #my @out_lines;

    my @well_keys = sort keys %$clip;

    return;
}

sub create_output_file {
    my $file_name = shift;
    my $lines = shift;

    my $tab_filename = ($ENV{'LIMS2_TEMP'} // '/var/tmp') . '/' . $file_name;
    open my $tab_fh, '>', $tab_filename
        or die "Can't open file $tab_filename: $! \n";
    print $tab_fh $lines;
    close $tab_fh
        or croak "Unable to close '$tab_filename' after writing: $!";

    return;
}


sub generate_crispr_output {
    my $out_rows = shift;
    my $headers = generate_crispr_headers();

    my $op;
    $op = "Sequencing primers\n";
    $op .= "WARNING: These primers are for sequencing a PCR product - no genomic check has been applied to these primers\n\n";
    $op .= $$headers . "\n";
    foreach my $row ( sort @{$out_rows} ) {
         $op .= $row . "\n";
    }
    $op .= "End of File\n";
    return $op;
}

sub generate_single_crispr_output {
    my $out_rows = shift;
    my $headers = generate_single_crispr_headers();

    my $op;
    $op = "Sequencing primers\n";
    $op .= "WARNING: These primers are for sequencing a PCR product - no genomic check has been applied to these primers\n\n";
    $op .= $$headers . "\n";
    foreach my $row ( sort @{$out_rows} ) {
         $op .= $row . "\n";
    }
    $op .= "End of File\n";
    return $op;
}

sub generate_pcr_output {
    my $out_rows = shift;
    my $headers = generate_pcr_headers();

    my $op;
    $op = "PCR crispr region primers\n";
    $op .= "Genomic specificity check has been applied to these primers\n\n";
    $op .= $$headers . "\n";
    foreach my $row ( sort @{$out_rows} ) {
         $op .= $row . "\n";
    }
    $op .= "End of File\n";
    return $op;
}

sub generate_genotyping_output {
    my $out_rows = shift;

    my $headers = generate_genotyping_headers();

    my $op;
    $op = "Genotyping primers\n";
    $op .= "Genomic specificity check has been applied to these primers\n\n";
    $op .= $$headers . "\n";
    foreach my $row ( sort @{$out_rows} ) {
         $op .= $row . "\n";
    }
    $op .= "End of File\n";
    return $op;
}

sub generate_gene_symbols_cache {
    my $params = shift;

    my $wells = $params->{'wells'};
    my $design_data_cache = $params->{'design_data_cache'};
    my $model = $params->{'model'};
    my $species = $params->{'species'};

    my $gene_cache;
    foreach my $well ( @{$wells} ) {
        my $well_id = $well->id;
        my $gene_id = $design_data_cache->{$well_id}->{'gene_id'};
        my $gene_symbol;


        if ( $gene_id ) {
            if ( $gene_cache->{$gene_id} ) {
                $gene_symbol = $gene_cache->{$gene_id};
            }
            else {
                my $gene_name_hr = $model->find_gene( {search_term => $gene_id, species => $species});
                if ( $gene_name_hr ) {
                    $gene_symbol = $gene_name_hr->{'gene_symbol'};
                }
                else {
                    $gene_symbol = '-';
                }
                $gene_cache->{$gene_id} = $gene_symbol;
            }
        }
        else {
            $gene_symbol = '-';
        }
        $design_data_cache->{$well_id}->{'gene_symbol'} = $gene_symbol;
    }
    return $design_data_cache;
}


sub generate_single_crispr_headers {

    my @crispr_headings = (qw/
        well_name
        gene_symbol
        design_id
        strand
        crispr_id
        SF1
        SF1_strand
        SF1_length
        SF1_gc_content
        SF1_tm
        SR1
        SR1_strand
        SF1_length
        SR1_gc_content
        SR1_tm
        crispr_id
        crispr_seq
    /);

    my $headers = join ',', @crispr_headings;

    return \$headers;
}

sub generate_crispr_headers {

    my @crispr_headings = (qw/
        well_name
        gene_symbol
        design_id
        strand
        crispr_pair_id
        SF1
        SF1_strand
        SF1_length
        SF1_gc_content
        SF1_tm
        SR1
        SR1_strand
        SF1_length
        SR1_gc_content
        SR1_tm
        crispr_left_id
        crispr_left_seq
        crispr_right_id
        crispr_right_seq
    /);

    my $headers = join ',', @crispr_headings;

    return \$headers;
}

sub generate_pcr_headers {

    my @pcr_headings = (qw/
        well_name
        gene_symbol
        design_id
        strand
        PF1
        PF1_strand
        PF1_length
        PF1_gc_content
        PF1_tm
        PR1
        PR1_strand
        PF1_length
        PR1_gc_content
        PR1_tm
        PF2
        PF2_strand
        PF2_length
        PF2_gc_content
        PF2_tm
        PR2
        PR2_strand
        PF2_length
        PR2_gc_content
        PR2_tm
    /);

    my $headers = join ',', @pcr_headings;

    return \$headers;
}

sub generate_genotyping_headers {
    my $design_oligos = shift;

    my @genotyping_headings = (qw/
        well_name
        gene_symbol
        design_id
        strand
        GF1
        GF1_strand
        GF1_length
        GF1_gc_content
        GF1_tm
        GR1
        GR1_strand
        GF1_length
        GR1_gc_content
        GR1_tm
        GF2
        GF2_strand
        GF2_length
        GF2_gc_content
        GF2_tm
        GR2
        GR2_strand
        GF2_length
        GR2_gc_content
        GR2_tm
        Design_A
        Oligo_A
        Design_B
        Oligo_B
    /);

    my $headers = join ',', @genotyping_headings;

    return \$headers;
}

sub generate_crispr_pair_id_cache {
    my $model = shift;
    my $plate_crispr_left = shift;
    my $plate_crispr_right = shift;

    my $crispr_left_rs;
    my $crispr_right_rs;
    my $cache = ();

    if ( $plate_crispr_left ) {
       $crispr_left_rs = $model->schema->resultset( 'Plate' )->search({
            'name' => $plate_crispr_left,
       });
    }
    if ( $plate_crispr_right ) {
        $crispr_right_rs = $model->schema->resultset( 'Plate' )->search({
            'name' => $plate_crispr_right,
            });
    }
    else {
        return $cache;
    }

    my $crispr_left_r = $crispr_left_rs->first;
    my $crispr_left_plate_name = $crispr_left_r->name;
    $logger->info( 'Crispr Left plate name retrieved: ' . $crispr_left_plate_name );

    my $crispr_right_r = $crispr_right_rs->first;
    my $crispr_right_plate_name = $crispr_right_r->name;
    $logger->info( 'Crispr Right plate name retrieved: ' . $crispr_right_plate_name );

    my @left_wells = $crispr_left_r->wells->all;
    $logger->info( 'Plate: ' . $crispr_left_plate_name .' - ' . @left_wells . ' wells retrieved');
    my @right_wells = $crispr_right_r->wells->all;
    $logger->info( 'Plate: ' . $crispr_right_plate_name .' - ' . @right_wells . ' wells retrieved');


    foreach my $left_well ( @left_wells ) {
        my $left_crispr_data;
        my $left_process_crispr = $left_well->process_output_wells->first->process->process_crispr;
        if ( $left_process_crispr ) {
            $left_crispr_data = $left_process_crispr->crispr->as_hash;
        }
        my $right_crispr_data;
        my ($right_well) = grep { $_->name eq $left_well->name } @right_wells;

        my $right_process_crispr = $right_well->process_output_wells->first->process->process_crispr;
        if ( $right_process_crispr ) {
            $right_crispr_data = $right_process_crispr->crispr->as_hash;
        }

        my $left_crispr_id = $left_crispr_data->{'id'};
        my $right_crispr_id = $right_crispr_data->{'id'};

        my $crispr_pair = $model->schema->resultset( 'CrisprPair' )->find({
           'left_crispr_id' => $left_crispr_id,
           'right_crispr_id' => $right_crispr_id,
        });


        $cache->{$left_well->name} = {
            'crispr_pair_id' =>  $crispr_pair->id,
        };
    }
    return $cache;

}
