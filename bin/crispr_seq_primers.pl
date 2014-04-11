#! /usr/bin/perl

use LIMS2::Model;
use feature qw/ say /;
use strict;
use warnings;
use Try::Tiny;

use LIMS2::Model::Util::OligoSelection;

=head filter_crisprs_coding

This file writes out genotyping_primers.csv

Calls the OligoSelection module of LIMS2 and formats the resulting hash
as a csv file.

=cut

my $model = LIMS2::Model->new( { user => 'webapp', audit_user => $ENV{USER} .'@sanger.ac.uk' } );

# Probably get the design ids from a file or the command line or directly from the database.
#my @design_id_list = ('1002436');
#my $crispr_pair_id =  19768;
my $species = 'Human';

my %primer_clip;
my $primer_results;
my $crispr_results;


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

my $design_row;
#while ($design_row = $crispr_design_rs->next ) {
foreach my $well ( @wells ) {
#    my $design_well = $well->design;
    my $well_id = $well->id;
    my $design_id = $design_data_cache->{$well_id}->{'design_id'};
    my $gene_id = $design_data_cache->{$well_id}->{'gene_id'};
    my $well_name = $well->name;
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
my ($crispr_design) = $model->schema->resultset('CrisprDesign')->search ( { 'design_id' => $design_id } );
    my $crispr_pair_id = $crispr_design->crispr_pair_id;
    say $design_id . ' - crispr_pair_id: ' . $crispr_pair_id; 

#    $primer_results = LIMS2::Model::Util::OligoSelection::pick_genotyping_primers( $model->schema, $design_id, $species )
    $primer_results = 'not yet implemented';
    $primer_clip{$well_name}{'genotyping'} = $primer_results;
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
foreach my $well_name ( keys %primer_clip ) {
    my @out_vals = (
        $well_name,
        $primer_clip{$well_name}{'gene_name'},
        $primer_clip{$well_name}{'design_id'},
        $primer_clip{$well_name}{'strand'},
        $primer_clip{$well_name}{'pair_id'});
    push (@out_vals, (
        $primer_clip{$well_name}->{'crispr_primers'}->{'left'}->{'left_0'}->{'seq'},
        $primer_clip{$well_name}->{'crispr_primers'}->{'left'}->{'left_0'}->{'location'}->{'_strand'},
        $primer_clip{$well_name}->{'crispr_primers'}->{'left'}->{'left_0'}->{'length'},
        $primer_clip{$well_name}->{'crispr_primers'}->{'left'}->{'left_0'}->{'gc_content'},
        $primer_clip{$well_name}->{'crispr_primers'}->{'left'}->{'left_0'}->{'melting_temp'},
    ));
    push (@out_vals, (
        $primer_clip{$well_name}->{'crispr_primers'}->{'right'}->{'right_0'}->{'seq'},
        $primer_clip{$well_name}->{'crispr_primers'}->{'right'}->{'right_0'}->{'location'}->{'_strand'},
        $primer_clip{$well_name}->{'crispr_primers'}->{'right'}->{'right_0'}->{'length'},
        $primer_clip{$well_name}->{'crispr_primers'}->{'right'}->{'right_0'}->{'gc_content'},
        $primer_clip{$well_name}->{'crispr_primers'}->{'right'}->{'right_0'}->{'melting_temp'},
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

my $headers = generate_headers();

my $tab_filename = 'crispr__primers.csv';
open( my $tab_fh, '>', $tab_filename )
    or die "Can't open file $tab_filename: $! \n";
print $tab_fh "Sequencing primers\n";
print $tab_fh "WARNING: These primers are for sequencing a PCR product - no genomic check has been applied to these primers\n\n";
print $tab_fh $$headers . "\n";
foreach my $row ( sort @out_rows ) {
    print $tab_fh $row;
    print $tab_fh "\n";
}
print $tab_fh "End of File\n";
close $tab_fh;

say 'Done' . "\n";

##
#

sub generate_headers {

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


=head
    # Find a left high ranking primer and select the corresponding right primer
    foreach my $exon_rank ( sort keys %{$clip->{$gene_symbol}} ) {
        if ( $clip->{$gene_symbol}->{$exon_rank}->{'5_prime'} ) {
            print $tab_fh "\n";
            foreach my $line ( @{$clip->{$gene_symbol}->{$exon_rank}->{'5_prime'}} ) {
                if ( defined($line) and length($line) ) {
                    print $tab_fh $line . "\n";
                }
                else
                {
                    print $tab_fh 'No results with selected parameters' . "\n";
                }
            }
        }
        else
        {
            print $tab_fh 'No results with selected parameters for this exon' . "\n";
        }
    }
    print $tab_fh "\n";
    print $tab_fh "$gene_symbol: 3_prime\n"; 
    print $tab_fh $headers_3_tsv . "\n";
    
    foreach my $exon_rank ( sort keys %{$clip->{$gene_symbol}} ) {
        if ( $clip->{$gene_symbol}->{$exon_rank}->{'3_prime'} ) {
            print $tab_fh "\n";
            foreach my $line ( @{$clip->{$gene_symbol}->{$exon_rank}->{'3_prime'}} ) {
                if ( defined($line) and length($line) ) {
                print $tab_fh $line . "\n"
                }
                else
                {
                    print $tab_fh 'No results with selected parameters' . "\n";
                }
            }
        }
        else
        {
            print $tab_fh 'No results with selected parameters for this exon' . "\n";
        }
}
=cut

