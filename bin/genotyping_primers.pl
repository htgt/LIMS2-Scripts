#! /usr/bin/perl

use LIMS2::Model;
use feature qw/ say /;
use strict;
use warnings;
use Try::Tiny;
use Carp;

use LIMS2::Model::Util::OligoSelection;

=head filter_crisprs_coding

This file writes out genotyping_primers.csv

Calls the OligoSelection module of LIMS2 and formats the resulting hash
as a csv file.

=cut

my $model = LIMS2::Model->new( { user => 'webapp', audit_user => $ENV{USER} .'@sanger.ac.uk' } );

# Probably get the design ids from a file or the command line or directly from the database.
my $species = 'Human';

my $plate_rs = $model->schema->resultset( 'Plate' )->search(
    {
        'name' => 'HG1'
    },
);


my $plate = $plate_rs->first;
my $plate_name = $plate->name;
print 'Plate selected: ' . $plate_name . "\n";

my @wells = $plate->wells->all;

my $well_count = @wells;
say 'Processing crispr primers for ' . $well_count . ' wells:';
my @well_id_list;

foreach my $well ( @wells ) {
    push @well_id_list, $well->id;
}

my $design_data_cache = $model->create_design_data_cache( \@well_id_list );

say 'Created design well data cache';


my $gene_cache;

run_primers();

say 'Done' . "\n";

##
#
sub run_primers {
    my $lines;
    say 'Generating gene symbol cache';

    generate_gene_symbols_cache();
    say 'Preparing crispr primers';
    my $out_rows = prepare_crispr_primers();
    say 'Generating crispr primer output file';
    $lines = generate_crispr_output( $out_rows );
    create_output_file( 'crispr_primers.csv', $lines );
    say 'Generating genotyping primers';
    $out_rows = prepare_genotyping_primers();
    $lines = generate_genotyping_output( $out_rows );
    say 'Generating genotyping primer output file';
    create_output_file( 'genotpying_primers.csv' ,$lines );

    return;
}


sub prepare_genotyping_primers {

    my %primer_clip;

    my $design_row;
    foreach my $well ( @wells ) {
        my $well_id = $well->id;
        my $design_id = $design_data_cache->{$well_id}->{'design_id'};
        my $gene_id = $design_data_cache->{$well_id}->{'gene_id'};
        my $well_name = $well->name;
        my $gene_name;

        $gene_name = $design_data_cache->{$well_id}->{'gene_symbol'};

        say $design_id . "\t(" . $gene_name . ')';

        my ($genotyping_primers, $genotyping_mapped, $chr_strand, $design_oligos)
            = LIMS2::Model::Util::OligoSelection::pick_genotyping_primers( {
                schema => $model->schema,
                design_id => $design_id,
                well_id => $well,
                species => $species,
            });
        $primer_clip{$well_name}{'gene_name'} = $gene_name;
        $primer_clip{$well_name}{'design_id'} = $design_id;
        $primer_clip{$well_name}{'strand'} = $chr_strand;

        $primer_clip{$well_name}{'genotyping_primers'} = $genotyping_mapped;
        $primer_clip{$well_name}{'design_oligos'} = $design_oligos;
    }
    my @out_rows;
    my $primer_type = 'genotyping_primers';
    foreach my $well_name ( keys %primer_clip ) {
        my @out_vals = (
            $well_name,
            $primer_clip{$well_name}{'gene_name'},
            $primer_clip{$well_name}{'design_id'},
            $primer_clip{$well_name}{'strand'},
        );
        # Take the two highest ranking primers
        my ($rank_a, $rank_b) = get_best_two_primer_ranks( $primer_clip{$well_name}->{$primer_type}->{'left'} );
        # only need left or right as both will be the same for rank purposes
        # GF1/GR1
        if ($primer_clip{$well_name}{'strand'} eq 'plus') {
            push (@out_vals, 
                data_to_push(\%primer_clip, $primer_type, $well_name, 'left', $rank_a // '99')
            );
            push (@out_vals, (
                data_to_push(\%primer_clip, $primer_type, $well_name, 'right', $rank_a // '99')
            ));
            # GF2/GR2
            push (@out_vals, 
                data_to_push(\%primer_clip, $primer_type, $well_name, 'left', $rank_b // '99')
            );
            push (@out_vals, (
                data_to_push(\%primer_clip, $primer_type, $well_name, 'right', $rank_b // '99')
            ));
        }
        elsif ($primer_clip{$well_name}{'strand'} eq 'minus') {
            push (@out_vals, 
                data_to_push(\%primer_clip, $primer_type, $well_name, 'right', $rank_a // '99')
            );
            push (@out_vals, (
                data_to_push(\%primer_clip, $primer_type, $well_name, 'left', $rank_a // '99')
            ));
            # GF2/GR2
            push (@out_vals, 
                data_to_push(\%primer_clip, $primer_type, $well_name, 'right', $rank_b // '99')
            );
            push (@out_vals, (
                data_to_push(\%primer_clip, $primer_type, $well_name, 'left', $rank_b // '99')
            ));

        }
        else {
            say 'Chromosome strand not defined.' ;
        }
        push (@out_vals,
            $primer_clip{$well_name}->{'design_oligos'}->{'5F'}->{'seq'},
            $primer_clip{$well_name}->{'design_oligos'}->{'3R'}->{'seq'},
        );

        my $csv_row = join( ',' , @out_vals);
        push @out_rows, $csv_row;
    }

    return \@out_rows;
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

sub prepare_crispr_primers {
    my %primer_clip;

    my $design_row;
    foreach my $well ( @wells ) {
        my $well_id = $well->id;
        my $design_id = $design_data_cache->{$well_id}->{'design_id'};
        my $gene_id = $design_data_cache->{$well_id}->{'gene_id'};
        my $well_name = $well->name;
        my $gene_name;

        $gene_name = $design_data_cache->{$well_id}->{'gene_symbol'};

        my ($crispr_design) = $model->schema->resultset('CrisprDesign')->search ( { 'design_id' => $design_id } );
        my $crispr_pair_id = $crispr_design->crispr_pair_id;
        say "$design_id\t$gene_name\tcrispr_pair_id:\t$crispr_pair_id";

        my ($crispr_results, $crispr_primers, $chr_strand) = LIMS2::Model::Util::OligoSelection::pick_crispr_primers( {
                schema => $model->schema,
                design_id => $design_id,
                crispr_pair_id => $crispr_pair_id,
                species => $species,
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

    return \@out_rows;
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

    foreach my $well ( @wells ) {
        my $well_id = $well->id;
        my $gene_id = $design_data_cache->{$well_id}->{'gene_id'};
        my $gene_name;


        if ( $gene_id ) {
            if ( $gene_cache->{$gene_id} ) {
                $gene_name = $gene_cache->{$gene_id};
            }
            else {
                $gene_name = $model->get_gene_symbol_for_gene_id( $gene_id, $species);
                $gene_cache->{$gene_id} = $gene_name;
            }
        }
        else {
            $gene_name = '-';
        }
        $design_data_cache->{$well_id}->{'gene_symbol'} = $gene_name;
    }
    return;
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

sub generate_genotyping_headers {

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
        5F
        3R
    /);

    my $headers = join ',', @genotyping_headings;

    return \$headers;
}
