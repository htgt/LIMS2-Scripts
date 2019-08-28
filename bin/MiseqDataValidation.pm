package MiseqDataValidation;

use strict;
use warnings;

use Sub::Exporter -setup => {
    exports => [qw(miseq_data_validator
                    get_all_miseq_plates
                    get_miseq_exps 
                    find_lost_wells 
                    get_well_ids 
                    get_well_name 
                    check_for_lost_well 
                    construct_file_path 
                    create_report_file 
                    write_report 
                )]
};

use LIMS2::Model;
use LIMS2::Model::Util::Miseq qw(convert_index_to_well_name);

my $schema = LIMS2::Model->new( user => 'lims2' )->schema;
my $plate_rs = $schema->resultset('Plate');
my $miseq_plate_rs = $schema->resultset('MiseqPlate');
my $miseq_exp_rs = $schema->resultset('MiseqExperiment');
my $miseq_well_exp_rs = $schema->resultset('MiseqWellExperiment');
my $well_rs = $schema->resultset('Well');
my $well_not_found = 1;

sub miseq_data_validator {
    #my @plates = get_all_miseq_plates();
    my $plate = $plate_rs->find({'name' => 'Miseq_010'})->as_hash;
    #foreach my $plate (@plates) {
        my %mismatch_cases = ();
        print "Scanning $plate->{name}\n";
        my @exps = get_miseq_exps($plate->{id});
        foreach my $exp (@exps) {
            if (my $lost_wells = find_lost_wells($exp->{id}, $plate->{name}, $exp->{name})) {
                $mismatch_cases{$exp->{name}} = $lost_wells;
            }
        }
        if (keys %mismatch_cases) {
            create_report_file($plate->{name}, \%mismatch_cases);
        }
    #}
}

sub get_all_miseq_plates {
    my @all_miseq_plates = map {$_->as_hash} $plate_rs->search({ type_id => 'MISEQ'});
    return @all_miseq_plates;
}

sub get_miseq_exps {
    my $plate_id = shift;
    #my $miseq_plate_id = $miseq_plate_rs->find({ plate_id => $plate_id })->id;
    my $miseq_plate_id = $miseq_plate_rs->find({ plate_id => $plate_id }, {column => 'id'})->id;
    my @miseq_exps = map {$_->as_hash} $miseq_exp_rs->search({ miseq_id => $miseq_plate_id });
    return @miseq_exps;
}

sub find_lost_wells {
    my ($exp_id, $plate_name, $exp_name) = @_;
    my @lost_wells = ();
    my @well_ids = get_well_ids($exp_id);
    my %existing_well_names = map {get_well_name($_) => 1} @well_ids;
    foreach (my $i = 1; $i < 385; $i++) {
        if (my $lost_well = check_for_lost_well($plate_name, $i, $exp_name, \%existing_well_names)) {
            push @lost_wells, $lost_well;
        }
    }
    return \@lost_wells;
}

sub get_well_ids {
    my $exp_id = shift;
    #my @well_ids = map {$_->well_id} $miseq_well_exp_rs->search({ miseq_exp_id => $exp_id});
    my @well_ids = map {$_->well_id} $miseq_well_exp_rs->search({ miseq_exp_id => $exp_id}, {column => 'well_id'});
    return @well_ids;
}

sub get_well_name {
    my $well_id = shift;
    #my $well_name = $well_rs->find({ id => $well_id })->name;
    my $well_name = $well_rs->find({ id => $well_id }, {column => 'name'})->name;
    return $well_name || $well_not_found;
}

sub check_for_lost_well {
    my ($plate_name, $i, $exp_name, $existing_well_names) = @_;
    my $file_path = construct_file_path($plate_name, $i, $exp_name);
    my $well_name_to_check = convert_index_to_well_name($i);
    if (-e $file_path) {
        if (! exists($existing_well_names->{$well_name_to_check})) {
            return $well_name_to_check;
        }
    }
    return;
}

sub construct_file_path {
    my ($plate_name, $i, $exp_name) = @_;
    my $file_path = "/warehouse/team229_wh01/lims2_managed_miseq_data/$plate_name/S${i}_exp$exp_name/CRISPResso_on_${i}_S${i}_L001_R1_001_${i}_S${i}_L001_R2_001/Quantification_of_editing_frequency.txt";
    return $file_path;
}

sub create_report_file {
    my ($plate_name, $mismatch_cases) = @_;
    my $report_filename = "${plate_name}_lost_wells.txt";
    open(my $report_fh, '>', $report_filename) or die "Couldn't open $report_filename: $!";
    write_report($report_fh, $mismatch_cases);
    close $report_fh;
    return;
}

sub write_report {
    my ($report_fh, $mismatch_cases) = @_;
    while (my ($exp, $lost_wells) = each %$mismatch_cases) {
        print $report_fh "$exp: @$lost_wells\n";
    }
    return;
}

1;
