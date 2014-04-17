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
my $species_param = 'Human'; # default is Human TODO: assembly?
my $crispr_type = 'pair';
my $format_param = '';
my @repeat_mask_param;

GetOptions(
    'plate=s'       => \$plate_name_param,
    'well=s'        => \$plate_well_param,
    'repeat_mask=s' => \@repeat_mask_param,
    'species=s'     => \$species_param,
    'crispr_type=s' => \$crispr_type,
    'format=s'      => \$format_param,
)
or die usage_message();;

if ($plate_name_param eq '') {
    die usage_message()
}

if ( scalar @repeat_mask_param == 0 ) {
    push  @repeat_mask_param,'NONE';
}

my $data_clip;
if ( $plate_well_param  eq '') {
    if ( $crispr_type eq 'single' ) {
        $data_clip = primers_for_single_crispr_plate( $plate_name_param, \@repeat_mask_param, $species_param );

    }
    else { 
        $data_clip = primers_for_plate( $plate_name_param, \@repeat_mask_param , $species_param );
    }
}
else
{
    if ( $crispr_type eq 'single' ) {
        # $plate_well_param is optional here
        $data_clip = primers_for_single_crispr_plate( $plate_name_param, \@repeat_mask_param, $species_param, $plate_well_param );
    }
    else {
        $data_clip = primers_for_plate_well( $plate_name_param, $plate_well_param, \@repeat_mask_param, $species_param );
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
    [--well=well_name]
    [--crispr_type=(single | pair)]
    [--species=(Human | Mouse)]
    [--repeat_mask=TRF [--repeat_mask=...] (default: NONE)]
    [--format=fsa]

Optional parameters in square brackets
Default species is Human
Default crispr_type is pair

END_DIE
}

# subroutines start here
#
sub primers_for_single_crispr_plate {
    my $plate_name_input = shift;
    my $repeat_mask = shift;
    my $species_input = shift;
    my $plate_well_input = shift;

    $logger->info("Starting primer generation for plate $plate_name_input");
    $logger->info("Primers for sequencing and pcr of a single crispr per well" );
    my $rpt_string = join( ',', @$repeat_mask);
    $logger->info("Using sequence repeat mask of: $rpt_string");

    my $model = LIMS2::Model->new( { user => 'webapp', audit_user => $ENV{USER} .'@sanger.ac.uk' } );

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


    my $sc_params = {
            'wells' => \@wells,
            'design_data_cache' => $design_data_cache,
            'model' => $model,
            'species' => $species,
            'plate_name' => $plate_name,
            'repeat_mask' => $repeat_mask,
        };

    if ( $plate_well_input ) {
        $sc_params->{'plate_well'} = $plate_well_input;
    }

    my $clip = run_single_crispr_primers( $sc_params );

    $logger->info( 'Done' );
    return $clip;
}

sub primers_for_plate {
    my $plate_name_input = shift;
    my $repeat_mask = shift;
    my $species_input = shift;

    $logger->info("Starting primer generation for plate $plate_name_input");
    my $rpt_string = join( ',', @$repeat_mask);
    $logger->info("Using sequence repeat mask of: $rpt_string");

    my $model = LIMS2::Model->new( { user => 'webapp', audit_user => $ENV{USER} .'@sanger.ac.uk' } );

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

    foreach my $well ( @wells ) {
        push @well_id_list, $well->id;
    }

    my $design_data_cache = $model->create_design_data_cache(
            \@well_id_list,
        );

    $logger->debug( 'Created design well data cache' );



    run_primers({
            'wells' => \@wells,
            'design_data_cache' => $design_data_cache,
            'model' => $model,
            'species' => $species,
            'plate_name' => $plate_name,
            'repeat_mask' => $repeat_mask,
        });

    $logger->info( 'Done' );
    return;
}

sub primers_for_plate_well {
    my $plate_name_input = shift;
    my $plate_well_input = shift;
    my $repeat_mask=shift;
    my $species_input = shift;

    $logger->info("Starting primer generation for plate $plate_name_input well $plate_well_input");
    my $rpt_string = join( ',', @$repeat_mask);
    $logger->info("Using sequence repeat mask of: $rpt_string");

    my $model = LIMS2::Model->new( { user => 'webapp', audit_user => $ENV{USER} .'@sanger.ac.uk' } );

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
    $logger->info( 'Processing crispr primers for 1 well of ' . $well_count . ' wells:');
    my @well_id_list;

    foreach my $well ( @wells ) {
        if ( $well->name eq $plate_well_input ) {
            push @well_id_list, $well->id;
        }
    }

    @wells = grep { $_->name =~ /$plate_well_input/  } @wells;

    my $design_data_cache = $model->create_design_data_cache(
            \@well_id_list,
        );

    $logger->debug( 'Created design well data cache' );



    run_primers({
            'wells' => \@wells,
            'design_data_cache' => $design_data_cache,
            'model' => $model,
            'species' => $species,
            'plate_name' => $plate_name,
            'repeat_mask' => $repeat_mask,
            'plate_well' => $plate_well_input,
        });

    $logger->info( 'Done' );
    return;
}

##
#

sub run_single_crispr_primers {
    my $params = shift;
    my $wells = $params->{'wells'};
    my $design_data_cache = $params->{'design_data_cache'};
    my $model = $params->{'model'};
    my $species = $params->{'species'};
    my $plate_name = $params->{'plate_name'};
    my $repeat_mask= $params->{'repeat_mask'};

    my $formatted_well = '_' . $params->{'plate_well'} // '';

    my $lines;
    $logger->debug( 'Generating gene symbol cache' );

    $design_data_cache = generate_gene_symbols_cache({
            'wells' => $wells,
            'design_data_cache' => $design_data_cache,
            'species' => $species,
            'model' => $model,
        });
    $logger->debug( 'Preparing crispr primers');
    my ($out_rows, $crispr_clip ) = prepare_single_crispr_primers({
            'model' => $model,
            'wells' => $wells,
            'design_data_cache' => $design_data_cache,
            'species' => $species,
            'repeat_mask' => $repeat_mask,
        });
    $logger->debug( 'Generating single crispr primer output file' );
    $lines = generate_single_crispr_output( $out_rows );
    create_output_file( $plate_name . $formatted_well . '_single_crispr_primers.csv', $lines );

    $logger->info( 'Generating PCR primers for crispr region' );
    ($out_rows, $crispr_clip) = prepare_pcr_primers({
            'model' => $model,
            'crispr_clip' => $crispr_clip,
            'wells' => $wells,
            'design_data_cache' => $design_data_cache,
            'species' => $species,
            'repeat_mask' => $repeat_mask,
        });
    $logger->debug( 'Generating pcr primer output file' );
    $lines = generate_pcr_output( $out_rows );
    create_output_file( $plate_name . $formatted_well . '_pcr_primers.csv', $lines );
    
    # now need to add back in genotyping primers in case short range primers not found
    # could just do this for the wells that failed in pcr long range.
    # For now just operate on all wells.
    $logger->info( 'Generating genotyping primers' );
    my $design_oligos;
    ($out_rows, $crispr_clip) = prepare_genotyping_primers({
            'model' => $model,
            'crispr_clip' => $crispr_clip,
            'wells' => $wells,
            'design_data_cache' => $design_data_cache,
            'species' => $species,
            'repeat_mask' => $repeat_mask,
        });

    $lines = generate_genotyping_output( $out_rows );
    $logger->info( 'Generating genotyping primer output file' );
    create_output_file( $plate_name . $formatted_well . '_genotyping_primers.csv' ,$lines );

    return $crispr_clip;
}

sub run_primers {
    my $params = shift;

    my $wells = $params->{'wells'};
    my $design_data_cache = $params->{'design_data_cache'};
    my $model = $params->{'model'};
    my $species = $params->{'species'};
    my $plate_name = $params->{'plate_name'};
    my $repeat_mask= $params->{'repeat_mask'};

    my $formatted_well = '_' . $params->{'plate_well'} // '';

    my $lines;
    $logger->debug( 'Generating gene symbol cache' );

    $design_data_cache = generate_gene_symbols_cache({
            'wells' => $wells,
            'design_data_cache' => $design_data_cache,
            'species' => $species,
            'model' => $model,
        });
    $logger->debug( 'Preparing crispr primers');
    my ($out_rows, $crispr_clip ) = prepare_crispr_primers({
            'model' => $model,
            'wells' => $wells,
            'design_data_cache' => $design_data_cache,
            'species' => $species,
            'repeat_mask' => $repeat_mask,
        });
    $logger->debug( 'Generating crispr primer output file' );
    $lines = generate_crispr_output( $out_rows );
    create_output_file( $plate_name . $formatted_well . '_crispr_primers.csv', $lines );

    $logger->info( 'Generating PCR primers for crispr region' );
    ($out_rows, $crispr_clip) = prepare_pcr_primers({
            'model' => $model,
            'crispr_clip' => $crispr_clip,
            'wells' => $wells,
            'design_data_cache' => $design_data_cache,
            'species' => $species,
            'repeat_mask' => $repeat_mask,
        });
    $logger->debug( 'Generating pcr primer output file' );
    $lines = generate_pcr_output( $out_rows );
    create_output_file( $plate_name . $formatted_well . '_pcr_primers.csv', $lines );


    $logger->info( 'Generating genotyping primers' );
    my $design_oligos;
    ($out_rows, $crispr_clip) = prepare_genotyping_primers({
            'model' => $model,
            'wells' => $wells,
            'design_data_cache' => $design_data_cache,
            'species' => $species,
            'repeat_mask' => $repeat_mask,
        });

    $lines = generate_genotyping_output( $out_rows );
    $logger->info( 'Generating genotyping primer output file' );
    create_output_file( $plate_name . $formatted_well . '_genotyping_primers.csv' ,$lines );

    return;
}


sub prepare_pcr_primers {
    my $params = shift;

    my $model = $params->{'model'};
    my $crispr_clip = $params->{'crispr_clip'};
    my $wells = $params->{'wells'};
    my $design_data_cache = $params->{'design_data_cache'};
    my $species = $params->{'species'};
    my $repeat_mask = $params->{'repeat_mask'};

    foreach my $well ( @{$wells} ) {
        my $well_id = $well->id;
        my $well_name = $well->name;

        my $gene_name = $design_data_cache->{$well_id}->{'gene_symbol'};

        $logger->info( 'crispr_pcr: ' . $well_name . "\t(" . $gene_name . ')' );

        my ($crispr_pcr_primers, $crispr_pcr_mapped)
            = LIMS2::Model::Util::OligoSelection::pick_crispr_PCR_primers( {
                schema => $model->schema,
                well_id => $well,
                crispr_primers => $crispr_clip->{$well_name},
                species => $species,
                repeat_mask => $repeat_mask,
            });
        $crispr_clip->{$well_name}->{'crispr_pcr_primers'} = $crispr_pcr_mapped;

    }
    my @out_rows;
    my $primer_type = 'crispr_pcr_primers';
    foreach my $well_name ( keys %{$crispr_clip} ) {
        my @out_vals = (
            $well_name,
            $crispr_clip->{$well_name}->{'gene_name'} // 'not set',
            $crispr_clip->{$well_name}->{'design_id'} // 'not set',
            $crispr_clip->{$well_name}->{'strand'} // 'not set',
        );
        # Take the two highest ranking primers
        my ($rank_a, $rank_b) = get_best_two_primer_ranks( $crispr_clip->{$well_name}->{$primer_type}->{'left'} );
        # only need left or right as both will be the same for rank purposes
        if ($crispr_clip->{$well_name}->{'strand'} eq 'plus') {
            # first set
            push (@out_vals,
                data_to_push($crispr_clip, $primer_type, $well_name, 'left', $rank_a // '99')
            );
            push (@out_vals, (
                data_to_push($crispr_clip, $primer_type, $well_name, 'right', $rank_a // '99')
            ));
            # second set
            push (@out_vals,
                data_to_push($crispr_clip, $primer_type, $well_name, 'left', $rank_b // '99')
            );
            push (@out_vals, (
                data_to_push($crispr_clip, $primer_type, $well_name, 'right', $rank_b // '99')
            ));
        }
        elsif ($crispr_clip->{$well_name}->{'strand'} eq 'minus') {
            push (@out_vals,
                data_to_push($crispr_clip, $primer_type, $well_name, 'right', $rank_a // '99')
            );
            push (@out_vals, (
                data_to_push($crispr_clip, $primer_type, $well_name, 'left', $rank_a // '99')
            ));
            # GF2/GR2
            push (@out_vals,
                data_to_push($crispr_clip, $primer_type, $well_name, 'right', $rank_b // '99')
            );
            push (@out_vals, (
                data_to_push($crispr_clip, $primer_type, $well_name, 'left', $rank_b // '99')
            ));

        }
        else {
            $logger->error( 'Chromosome strand not defined.' );
        }
        # Anything after here for context only

        my $csv_row = join( ',' , @out_vals);
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

        my ($genotyping_primers, $genotyping_mapped, $chr_strand, $design_oligos)
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
    }
    my @out_rows;
    my $primer_type = 'genotyping_primers';
    foreach my $well_name ( keys %$primer_clip ) {
        my @out_vals = (
            $well_name,
            $primer_clip->{$well_name}->{'gene_name'},
            $primer_clip->{$well_name}->{'design_id'},
            $primer_clip->{$well_name}->{'strand'},
        );
        # Take the two highest ranking primers
        my ($rank_a, $rank_b) = get_best_two_primer_ranks( $primer_clip->{$well_name}->{$primer_type}->{'left'} );
        # only need left or right as both will be the same for rank purposes
        # GF1/GR1
        if ($primer_clip->{$well_name}->{'strand'} eq 'plus') {
            push (@out_vals,
                data_to_push($primer_clip, $primer_type, $well_name, 'left', $rank_a // '99')
            );
            push (@out_vals, (
                data_to_push($primer_clip, $primer_type, $well_name, 'right', $rank_a // '99')
            ));
            # GF2/GR2
            push (@out_vals,
                data_to_push($primer_clip, $primer_type, $well_name, 'left', $rank_b // '99')
            );
            push (@out_vals, (
                data_to_push($primer_clip, $primer_type, $well_name, 'right', $rank_b // '99')
            ));
        }
        elsif ($primer_clip->{$well_name}->{'strand'} eq 'minus') {
            push (@out_vals,
                data_to_push($primer_clip, $primer_type, $well_name, 'right', $rank_a // '99')
            );
            push (@out_vals, (
                data_to_push($primer_clip, $primer_type, $well_name, 'left', $rank_a // '99')
            ));
            # GF2/GR2
            push (@out_vals,
                data_to_push($primer_clip, $primer_type, $well_name, 'right', $rank_b // '99')
            );
            push (@out_vals, (
                data_to_push($primer_clip, $primer_type, $well_name, 'left', $rank_b // '99')
            ));

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

        my $csv_row = join( ',' , @out_vals);
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
    my $repeat_mask = $params->{'repeat_mask'};

    my %primer_clip;

    my $design_row;
    foreach my $well ( @{$wells} ) {
        my $well_id = $well->id;
        my $design_id = $design_data_cache->{$well_id}->{'design_id'};
        my $gene_id = $design_data_cache->{$well_id}->{'gene_id'};
        my $well_name = $well->name;
        my $gene_name;

        $gene_name = $design_data_cache->{$well_id}->{'gene_symbol'};
$DB::single=1;

        my ($crispr_left, $crispr_right) = $well->left_and_right_crispr_wells;
        my $crispr_id = $crispr_left->crispr->id;
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
            data_to_push(\%primer_clip, $primer_type, $well_name, 'left', $rank)
        ));
        push (@out_vals, (
            data_to_push(\%primer_clip, $primer_type, $well_name, 'right', $rank)
        ));
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

    my %primer_clip;

    my $design_row;
    foreach my $well ( @{$wells} ) {
        my $well_id = $well->id;
        my $design_id = $design_data_cache->{$well_id}->{'design_id'};
        my $gene_id = $design_data_cache->{$well_id}->{'gene_id'};
        my $well_name = $well->name;
        my $gene_name;

        $gene_name = $design_data_cache->{$well_id}->{'gene_symbol'};

        my ($crispr_design) = $model->schema->resultset('CrisprDesign')->search ( { 'design_id' => $design_id } );
        my $crispr_pair_id = $crispr_design->crispr_pair_id;
        $logger->info( "$design_id\t$gene_name\tcrispr_pair_id:\t$crispr_pair_id" );

        my ($crispr_results, $crispr_primers, $chr_strand) = LIMS2::Model::Util::OligoSelection::pick_crispr_primers( {
                schema => $model->schema,
                design_id => $design_id,
                crispr_pair_id => $crispr_pair_id,
                species => $species,
                repeat_mask => $repeat_mask,
            });
        $primer_clip{$well_name}{'pair_id'} = $crispr_pair_id;
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
            $primer_clip{$well_name}{'pair_id'},
        );
        push (@out_vals, (
            data_to_push(\%primer_clip, $primer_type, $well_name, 'left', $rank)
        ));
        push (@out_vals, (
            data_to_push(\%primer_clip, $primer_type, $well_name, 'right', $rank)
        ));
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

=head generate_fsa_file
Generates data to output fasta format files for use with BlastN

=cut

sub generate_fsa_files {
    my $clip = shift;

    my @out_lines;

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
        my $gene_name;


        if ( $gene_id ) {
            if ( $gene_cache->{$gene_id} ) {
                $gene_name = $gene_cache->{$gene_id};
            }
            else {
                $gene_name = $model->find_gene( {search_term => $gene_id, species => $species});
                $gene_cache->{$gene_id} = $gene_name;
            }
        }
        else {
            $gene_name = '-';
        }
        $design_data_cache->{$well_id}->{'gene_symbol'} = $gene_name;
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
