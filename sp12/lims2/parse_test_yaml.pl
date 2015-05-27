#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Path::Class;
use YAML::Any qw(LoadFile);
use feature qw( say );
use Smart::Comments;

my $test_dir = dir( $ARGV[0] );

my $dir_name = $test_dir->absolute->stringify;

my @files = qw(
create_di_process.yaml
create_crispr_process.yaml
int_recom_process.yaml
2w_gateway_process.yaml
legacy_gateway_process.yaml
3w_gateway_process.yaml
cre_bac_recom_process.yaml
recombinase_process.yaml 
add_recombinase.yaml
final_pick_process.yaml
rearray_process.yaml
dna_prep_process.yaml
first_electroporation.yaml
second_electroporation.yaml
clone_pick_process.yaml
clone_pool_process.yaml
freeze_process.yaml
dist_qc_process.yaml
xep_pool_process.yaml
);

my %wells;
for my $file ( @files ) {
    my $data = LoadFile( $dir_name . '/' . $file  );

    for my $test ( keys %$data ) {
        my $datum = $data->{ $test };
        if ( my $output_wells = $datum->{output_wells} ) {
            for my $well ( @$output_wells ) {
                $wells{ $well->{plate_name} }{ $well->{well_name} }++;
            }
        }
        if ( my $output_wells = $datum->{input_wells} ) {
            for my $well ( @$output_wells ) {
                $wells{ $well->{plate_name} }{ $well->{well_name} }++;
            }
        }
    }
}

for my $plate ( keys %wells ) {
    say "$plate";
}

