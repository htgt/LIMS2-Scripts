#! /usr/bin/perl

# use LIMS2::Model;
# use feature qw/ say /;
use strict;
use warnings;

use LIMS2::Model;
use LIMS2::Model::Plugin::Project;
use LIMS2::ReportGenerator;
use Smart::Comments;
use TryCatch;
use feature qw( say );
use Pod::Usage;
use Getopt::Long;
use Log::Log4perl ':easy';



my $log_level = $WARN;
GetOptions(
    'help'                 => sub { pod2usage( -verbose => 1 ) },
    'man'                  => sub { pod2usage( -verbose => 2 ) },
    'debug'                => sub { $log_level = $DEBUG },
    'verbose'              => sub { $log_level = $INFO },
    'trace'                => sub { $log_level = $TRACE },
    'crispr_ep_file=s'     => \my $crispr_ep_file,
    'ep_pick_file=s'       => \my $ep_pick_file,
    'species=s'            => \my $species,
) or pod2usage(2);


my $model = LIMS2::Model->new( user => 'lims2' );



$species = 'Mouse' unless $species;
WARN( "Species: $species" );

$crispr_ep_file = "crispr_ep_report.csv" unless $crispr_ep_file;
WARN( "CRISPR_EP output file: $crispr_ep_file" );

$ep_pick_file = "ep_pick_report.csv" unless $ep_pick_file;
WARN( "EP_PICK output file: $ep_pick_file" );



my @crispr_ep_plates = $model->schema->resultset("Plate")->search(
    {
        species_id  => $species,
        type_id     => 'CRISPR_EP',
    }
);


open CRISPR_EP, ">", $crispr_ep_file or die $!;
open EP_PICK, ">", $ep_pick_file or die $!;


print CRISPR_EP join (',',
        'Plate ID', 'Plate Name',
        'Well ID', 'Well Name', 'Design ID', 'Design Type', 'Gene ID', 'Gene Symbol', 'Gene Sponsors', 'Genbank File',
        'Cassette', 'Cassette Resistance', 'Cassette Type', 'Backbone', 'Nuclease', 'Cell Line',
        'Picked Colonies', 'Remaining Stained Colonies', 'Remaining Unstained Colonies', 'Total Colonies',
        'Crispr Wells',
        'Created By','Created At', 'Report?', "\n"
    );

print EP_PICK join (',',
        'EP Plate ID', 'EP Well ID', 'EP Well Name','Plate ID', 'Plate Name',
        'Well ID', 'Well Name', 'Design ID', 'Design Type', 'Gene ID', 'Gene Symbol', 'Gene Sponsors',
        'Created By', 'Created At', 'Assay Pending', 'Assay Complete', 'Accepted?', 'Genbank File',
        'Cassette', 'Cassette Resistance', 'Recombinases', 'Cell Line',
        "Clone ID", "QC Pass", "Valid Primers", "QC Result URL", "Primer Bands", "\n"
    );


foreach my $crispr_ep (@crispr_ep_plates) {

    my $rs = $model->schema->resultset( 'PlateReport' )->search(
        {},
        {
            bind => [ $crispr_ep->id ],
        }
    );

    my @wells_data = @{ $rs->consolidate( $crispr_ep->id,
            [ 'well_qc_sequencing_result', 'well_colony_counts' ] ) };
    @wells_data = sort { $a->{well_name} cmp $b->{well_name} } @wells_data;


    while (@wells_data) {
        my $well_data = shift @wells_data;

        my @crispr_wells = map { $_->{plate_name} . '[' . $_->{well_name} . ']' }
            @{ $well_data->{crispr_wells}{crisprs} };

        my @data = (
            $crispr_ep->id,
            $crispr_ep->name,
            $well_data->{well_id},
            $well_data->{well_name},
            $well_data->{design_id},
            $well_data->{design_type},
            $well_data->{gene_ids},
            $well_data->{gene_symbols},
            $well_data->{sponsors},
            &well_eng_seq_link( $well_data ),
            $well_data->{cassette} // '',
            $well_data->{cassette_resistance} // '',
            $well_data->{cassette_promoter},
            $well_data->{backbone} // '',
            $well_data->{nuclease},
            $well_data->{cell_line},
            &colony_counts( $well_data->{well} ),
            '"' . join( ', ', @crispr_wells ) . '"',
            $well_data->{created_by},
            $well_data->{created_at},
            $well_data->{to_report} ? 'true' : 'false',
        );

        print CRISPR_EP join (',', @data, "\n");

    }



    my $children;

    for my $well ( $crispr_ep->wells ){
        foreach my $process ($well->child_processes){
            my $type = $process->type_id;
            next unless $type eq 'clone_pick';

            $children->{$type} ||= {};
            foreach my $output ($process->output_wells){
                my $plate = $output->plate;
                $children->{$type}->{$plate->id} = $plate;
            }
        }
    }


    foreach my $type (keys %{ $children || {} }){
        foreach my $child_plate_id (keys %{ $children->{$type}  }){

            my $ep_pick = $model->schema->resultset("Plate")->find(
                {
                    species_id  => $species,
                    id          => $child_plate_id,
                }
            );

            my $rs = $model->schema->resultset( 'PlateReport' )->search(
                {},
                {
                    bind => [ $ep_pick->id ],
                }
            );

            my @wells_data = @{ $rs->consolidate( $ep_pick->id,
                    [ 'well_qc_sequencing_result', 'well_colony_counts' ] ) };
            @wells_data = sort { $a->{well_name} cmp $b->{well_name} } @wells_data;


            while (@wells_data) {
                my $well_data = shift @wells_data;

                my @crispr_wells = map { $_->{plate_name} . '[' . $_->{well_name} . ']' }
                    @{ $well_data->{crispr_wells}{crisprs} };

                my @data = (
                    $crispr_ep->id,
                    $well_data->{well_ancestors}{CRISPR_EP}{well_id},
                    $well_data->{well_ancestors}{CRISPR_EP}{well_name},
                    $ep_pick->id,
                    $ep_pick->name,
                    $well_data->{well_id},
                    $well_data->{well_name},
                    $well_data->{design_id},
                    $well_data->{design_type},
                    $well_data->{gene_ids},
                    $well_data->{gene_symbols},
                    $well_data->{sponsors},

                    $well_data->{created_by},
                    $well_data->{created_at},
                    $well_data->{assay_pending},
                    $well_data->{assay_complete},
                    LIMS2::ReportGenerator->boolean_str( $well_data->{accepted} ),
                    &well_eng_seq_link( $well_data ),
                    $well_data->{cassette} // '',
                    $well_data->{cassette_resistance} // '',
                    $well_data->{recombinases} // '',
                    $well_data->{cell_line},

                    $ep_pick->name . '_' . $well_data->{well_name},
                    &well_qc_sequencing_result_data( $well_data->{well} ),
                    '"' .  &well_primer_bands_data( $well_data->{well} ) . '"',
                );

                print EP_PICK join (',', @data, "\n");

            }

        }
    }


}

close CRISPR_EP;
close EP_PICK;


sub colony_counts {
    my ( $well ) = @_;

    my @types =  map { $_->id } $model->schema->resultset('ColonyCountType')->search(
        {},
        { order_by => { -asc => 'id' } }
    );

    my $colony_counts = { map { $_ => '' } @types };

    foreach my $count ($well->well_colony_counts) {
        my $type = $count->colony_count_type_id;
        my $count = $count->colony_count;

        $colony_counts->{$type} = $count;

    }

    my $result = [
        $colony_counts->{'picked_colonies'},
        $colony_counts->{'remaining_stained_colonies'},
        $colony_counts->{'remaining_unstained_colonies'},
        $colony_counts->{'total_colonies'}
    ];

    return @{$result};
}


sub well_eng_seq_link {
    my ( $well_data ) = @_;

    # nonsense designs can not currently have eng seqs generated for them..
    return 'N/A' if $well_data->{design_type} eq 'nonsense';

    if ( $well_data->{well_id} ) {
        return 'http://www.sanger.ac.uk/htgt/lims2/public_reports/well_eng_seq/' . $well_data->{well_id};
    }

    return '-';
}


sub well_qc_sequencing_result_data {
    my ( $well ) = @_;

    my @qc_data = ( '', '', '' );
    if ( my $well_qc_sequencing_result = $well->well_qc_sequencing_result ) {

        @qc_data = (
            LIMS2::ReportGenerator->boolean_str( $well_qc_sequencing_result->pass ),
            '"' . $well_qc_sequencing_result->valid_primers . '"',
            $well_qc_sequencing_result->test_result_url,
        );
    }

    return @qc_data;
}


sub well_primer_bands_data {
    my ( $well ) = @_;

    my @primer_bands_data;
    for my $primer_band ( $well->well_primer_bands->all ) {
        push @primer_bands_data, $primer_band->primer_band_type_id . '(' . $primer_band->pass . ')';
    }

    return join( ', ', @primer_bands_data );
}


__END__

=head1 NAME

jogi_crispr_ep_report.pl - Make CRISPR_EP and EP_PICK well reports for Jogi ( Joachim Beig )

=head1 SYNOPSIS

  jogi_crispr_ep_report.pl [options]
      --help            Display a brief help message
      --man             Display the manual page
      --debug           Debug output
      --verbose         Verbose output
      --crispr_ep_file  File for the CRISPR_EP report ( default crispr_ep_report.csv )
      --ep_pick_file    File for the EP_PICK report ( default ep_pick_report.csv )
      --species         Species of targets ( default Mouse )


    This script will generate 2 reports:
     - a CRISPR_EP plate report
     - a EP_PICK plate report

    Originaly created by request of Jogi.

=head1 DESCRIPTION


=cut
