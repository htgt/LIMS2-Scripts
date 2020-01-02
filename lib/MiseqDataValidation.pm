package MiseqDataValidation;

use strict;
use warnings;

use Sub::Exporter -setup => {
    exports => [
        qw(
          miseq_data_validator
          get_all_miseq_plates
          check_plate
          get_mismatches
          get_miseq_exps
          find_lost_wells
          get_well_names
          check_for_lost_well
          construct_file_path
          create_report_file
          write_report
          )
    ]
};
use Net::OpenSSH;

use LIMS2::Model;
use LIMS2::Model::Util::Miseq qw(convert_index_to_well_name);

my $schema = LIMS2::Model->new( user => 'lims2' )->schema;

my $quant_file_root = '/warehouse/team229_wh01/lims2_managed_miseq_data';
my @quant_files;

sub miseq_data_validator {
    if ( $ENV{LIMS2_FILE_ACCESS_SERVER} ) {
        # This environment variable must contain a host for accessing remote files
        get_quant_files();
    }
    elsif ( !-e $quant_file_root ) {
        die "Quantification files won't be found";
    }
    my @plates = get_all_miseq_plates();
    foreach my $plate (@plates) {
        check_plate($plate);
    }
    return;
}

sub get_quant_files {
    print "Getting quantification files...\n";
    my $ssh = Net::OpenSSH->new( $ENV{LIMS2_FILE_ACCESS_SERVER} );
    @quant_files = $ssh->capture(
        "find $quant_file_root -name Quantification_of_editing_frequency.txt");
    chomp(@quant_files);
    return;
}

sub get_all_miseq_plates {
    my $plate_rs = $schema->resultset('Plate');
    my @all_miseq_plates =
      map { $_->as_hash } $plate_rs->search( { type_id => 'MISEQ' } );
    return @all_miseq_plates;
}

sub check_plate {
    my $plate = shift;
    print "Scanning $plate->{name}\n";
    my $plate_mismatches = get_mismatches( $plate->{id}, $plate->{name} );
    if ( keys %$plate_mismatches ) {
        create_report_file( $plate->{name}, $plate_mismatches );
    }
    return;
}

sub get_mismatches {
    my ( $plate_id, $plate_name ) = @_;
    my %plate_mismatches;
    my @miseq_exps = get_miseq_exps($plate_id);
    foreach my $miseq_exp (@miseq_exps) {
        if ( my @lost_wells = find_lost_wells( $miseq_exp, $plate_name ) ) {
            $plate_mismatches{ $miseq_exp->{name} } = \@lost_wells;
        }
    }
    return \%plate_mismatches;
}

sub get_miseq_exps {
    my $plate_id       = shift;
    my $miseq_plate_rs = $schema->resultset('MiseqPlate');
    my $miseq_exp_rs   = $schema->resultset('MiseqExperiment');
    if ( my $miseq_plate = $miseq_plate_rs->find( { plate_id => $plate_id } ) )
    {
        my $miseq_plate_id = $miseq_plate->id;
        my @miseq_exps     = map { $_->as_hash }
          $miseq_exp_rs->search( { miseq_id => $miseq_plate_id } );
        return @miseq_exps;
    }
    else {
        return;
    }
}

sub find_lost_wells {
    my ( $miseq_exp, $plate_name ) = @_;
    my @lost_wells;
    my $existing_well_names = get_well_names( $miseq_exp->{id} );
    foreach ( my $i = 1 ; $i < 385 ; $i++ ) {
        if (
            my $lost_well = check_for_lost_well(
                $plate_name, $i, $miseq_exp->{name}, $existing_well_names
            )
          )
        {
            push @lost_wells, $lost_well;
        }
    }
    return @lost_wells;
}

sub get_well_names {
    my $miseq_exp_id      = shift;
    my $miseq_well_exp_rs = $schema->resultset('MiseqWellExperiment');
    my @well_names        = map { $_->well->well_name }
      $miseq_well_exp_rs->search( { miseq_exp_id => $miseq_exp_id } );
    return \@well_names;
}

sub check_for_lost_well {
    my ( $plate_name, $i, $miseq_exp_name, $existing_well_names ) = @_;
    my $file_path = construct_file_path( $plate_name, $i, $miseq_exp_name );
    my $well_name_to_check = convert_index_to_well_name($i);
    if ( !grep { $well_name_to_check eq $_ } @$existing_well_names ) {
        if ( ( grep { $file_path eq $_ } @quant_files ) or ( -e $file_path ) ) {
            return $well_name_to_check;
        }
    }
    return;
}

sub construct_file_path {
    my ( $plate_name, $i, $miseq_exp_name ) = @_;
    my $file_path =
"$quant_file_root/$plate_name/S${i}_exp$miseq_exp_name/CRISPResso_on_${i}_S${i}_L001_R1_001_${i}_S${i}_L001_R2_001/Quantification_of_editing_frequency.txt";
    return $file_path;
}

sub create_report_file {
    my ( $plate_name, $plate_mismatches ) = @_;
    my $report_filename = "${plate_name}_lost_wells.txt";
    open( my $report_fh, '>', $report_filename )
      or die "Couldn't open $report_filename: $!";
    write_report( $report_fh, $plate_mismatches );
    close $report_fh or die "Couldn't close $report_filename: $!";
    return;
}

sub write_report {
    my ( $report_fh, $plate_mismatches ) = @_;
    while ( my ( $exp, $lost_wells ) = each %$plate_mismatches ) {
        print $report_fh "$exp: @$lost_wells\n";
    }
    return;
}

1;
