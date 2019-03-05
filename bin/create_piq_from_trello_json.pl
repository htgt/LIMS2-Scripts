#! /usr/bin/perl

use strict;
use warnings;

use LIMS2::Model;
use Try::Tiny;
use Getopt::Long;
use Data::Dumper;
use Moose;
use LIMS2::REST::Client;
use POSIX qw(strftime);
use feature qw(say);
use JSON;
use Text::CSV;

my $lims2_model = LIMS2::Model->new( user => 'lims2' );

sub generate_piq_well_names {
    my $start = shift;

    my @wells;
    my @letters = ('A'..'H');
    for (my $num = $start; $num <= $start + 2; $num++) {
        foreach my $letter (@letters) {
            my $name = $letter . sprintf("%02d", $num);
            push (@wells, $name);
        }
    }

    return @wells;
}

sub remove_well_leading_zero {
    my $name = shift;

    my @segs = split /_/, $name;
    my $plate = $segs[0];
    my $well = $segs[1] + 0;
    my $subclone = $segs[2];
    my $str;
    if ($subclone) {
        $str = $plate . '_' . $well . '_' . $subclone;
    } else {
        $str = $plate . '_' . $well;
    } 

    return $str;
}

sub update_miseq_experiments {
    my ($fp_name, $fp_id, $piq_name, $ep_info, @eps_piq) = @_;
    
    my ($plate_standard, $plate_ident, $clone) = $fp_name =~ /^(\D+)([\d]+)([_A-Z\d]+)$/;
    
    $plate_ident = sprintf("%04d", $plate_ident);
    my $epd_name = 'HUPEPD' . $plate_ident . $clone;
    my $msq_exp_rs = $lims2_model->schema->resultset('MiseqExperiment')->search({
        name => { -like => $epd_name . '%' },
    });
    
    if ($msq_exp_rs->count == 0) {
        print "\nNo Miseq Exp info for $epd_name - $fp_name";
    }
    
    while (my $msq_exp = $msq_exp_rs->next) {
        if ($ep_info->{$epd_name}) {
            my $exp_name = $msq_exp->name;
            #print "\nUpdating Miseq Exp: $exp_name - $ep_info->{$epd_name}->{experiment} : $ep_info->{$epd_name}->{gene}";
            $lims2_model->update_miseq_experiment({ 
                id              => $msq_exp->id,
                experiment_id   => $ep_info->{$epd_name}->{experiment},
                parent_plate_id => $fp_id,
                gene            => $ep_info->{$epd_name}->{gene},
            });
            
            my $piq_exp = {
                piq             => $piq_name,
                name            => $piq_name . '_' . $ep_info->{$epd_name}->{trivial},
                gene            => $ep_info->{$epd_name}->{gene},
                experiment_id   => $ep_info->{$epd_name}->{experiment},
            };
            push (@eps_piq, $piq_exp);
        } else {
            print "\nNo EP info for $epd_name - $fp_name";
        }
    }
    return @eps_piq;
}

my $dir;
my $dump;

GetOptions(
    'dir=s'  => \$dir,
    'dump'   => \$dump,
)
or die usage_message();;

opendir my $dh, $dir or die "Cannot open directory: $!";
my @paths = readdir $dh;
closedir $dh;

my $mapping;
my $mapping_file = $dir . 'PIQ_to_Miseq_Mapping.csv';
my $csv = Text::CSV->new ({ binary => 1, auto_diag => 1 });
open my $fh, $mapping_file or die "Cannot open mapping csv: $!";
$csv->column_names( $csv->getline($fh) );
while (my $row = $csv->getline_hr($fh)) {
    $mapping->{ $row->{piq_name} } = {
        miseq       => $row->{miseq_name},
        start_well  => $row->{miseq_start_well},
    };
}
close $fh;

my @files = grep (/.*\.json$/, @paths);
my $json = JSON->new;

my $name_regex = '^\ ?(\d+)\)\ ?([A-z0-9]+)\ +\S*\ *([A-Z][\dOI]{2}).*$';
# 2)HUPFP0033_2 RNF43_3 D06 WT  || 1)HUPFP0033_2 RNF43_3 D05 HET BN
#(2)(HUPFP0033_2)      (D06)    ||(1)(HUPFP0033_2)      (D05)

my @wells = generate_piq_well_names(1);

my $eps;
my @eps_piqs;
my $ep_file = $dir . 'Primary_EP.csv';
open my $epfh, $ep_file or die "Cannot open ep csv: $!";
$csv->column_names( $csv->getline($epfh) );
while (my $row = $csv->getline_hr($epfh)) {
    $eps->{ $row->{Cell_Number} } = {
        gene        => $row->{Gene_ID},
        experiment  => $row->{Experiment_ID},
        trivial     => $row->{Experiment},
    };
}
close $epfh;

my $unknown_entities;
my $piq_data;
my $miseq_data;

foreach my $json_filepath (@files) {
    my $miseq_check = 1;
    my $full_path = $dir . $json_filepath;
    my $json_text = do {
        open(my $json_fh, "<:encoding(UTF-8)", $full_path)
            or die("Can't open \$filename\": $!\n");
        local $/;
        my $json_data = <$json_fh>;
        close $json_fh;
        $json_data
    };

    my $card_data = $json->decode($json_text);

    my @name_split = split /\./, $json_filepath;
    my $piq_name = uc $name_split[0];
    my $piq_map = $mapping->{$piq_name};
    my $start_index = (split /A/, $piq_map->{start_well})[-1] + 0;

    my $miseq_plate_rs = $lims2_model->schema->resultset('Plate')->search({ name => $piq_map->{miseq} });
    if ($miseq_plate_rs->count == 0) {
        print "No Miseq plate record found for $piq_name : $piq_map->{miseq}";
        $miseq_check = 0;
    }
    my @miseq_ref = generate_piq_well_names($start_index);
    
    $piq_data->{$piq_name} = {
        name            => $piq_name,
        created_by      => 'pk8@sanger.ac.uk',
        species         => 'Human',
        created_at      => strftime("%Y-%m-%dT%H:%M:%S", localtime(time)),
        type            => 'PIQ',
    };
    
    my $arr_num = 0;
    if (@{$card_data->{checklists}[$arr_num]->{checkItems}} == 0) {
        $arr_num = 1;
    }


    my @parent_data = @{ $card_data->{checklists}[$arr_num]->{checkItems} };
    my @piq_wells;
    my @miseq_processes;
    foreach my $parent (@parent_data) {
        my $name_str = $parent->{name};
        my ($piq_inc, $plate_name, $well_name) = $name_str =~ m/$name_regex/;
        if ($plate_name) {
            $plate_name = remove_well_leading_zero($plate_name);
        } else { 
            push (@{$unknown_entities->{$piq_name}}, $name_str);
            next;
        }
        my $plate_rs = $lims2_model->schema->resultset('Plate')->search({ name => $plate_name });
        if ($plate_rs->count == 0) {
            $unknown_entities->{$plate_name} = "All - $piq_name: $name_str";
            next;
        }
        $well_name =~ s/O/0/g; 
        $well_name =~ s/I/1/g; 
        my $well_rs = $lims2_model->schema->resultset('Well')->search({ 
            name        => $well_name,
            plate_id    => { -in => $plate_rs->get_column('id')->as_query },
        });
        if ($well_rs->count == 0) {
            $unknown_entities->{$plate_name} = $well_name;
            next; 
        }
        
        my $piq_well_name = $wells[$piq_inc - 1];
        my $well_data = {
            well_name       => $piq_well_name,
            process_type    => 'dist_qc',
            parent_well     => $well_rs->first->name,
            parent_plate    => $plate_rs->first->name,
        };

        push (@piq_wells, $well_data);

        if ($miseq_check) { 
            @eps_piqs = update_miseq_experiments($plate_name, $plate_rs->first->id, $piq_name, $eps, @eps_piqs);

            my $miseq_well_name = $miseq_ref[$piq_inc - 1];

            my $process = {
                type => 'miseq_no_template',
                input_wells => [{
                    plate_name  => $piq_name,
                    well_name   => $piq_well_name,
                }],
                output_wells => [{
                    plate_name  => $piq_map->{miseq},
                    well_name   => $miseq_well_name,
                }],
            };
            push(@miseq_processes, $process);
        }
    }

    @piq_wells = sort { $a->{well_name} cmp $b->{well_name} } @piq_wells;
    $piq_data->{$piq_name}->{wells} = \@piq_wells;
    my $miseq_rs = $lims2_model->schema->resultset('Plate')->find({ name => $piq_map->{miseq} });
    $miseq_data->{$piq_name}->{$piq_map->{miseq}} = {
        piq         => $piq_name,
        miseq_dets  => $miseq_rs->miseq_details,
        wells       => \@miseq_processes,
    };

    my $egg = 1;
}

foreach my $piq_key (keys %{ $piq_data }) {
    my $piq_req = $piq_data->{$piq_key};
    my $curr_piq_rs = $lims2_model->create_plate($piq_req);
    print "\n~~~~~~~~~~~~~~~~~~~~~ Created PIQ plate: $piq_key ~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
}

foreach my $piq_key (keys %{ $miseq_data }) {
    foreach my $relation (values %{ $miseq_data->{$piq_key} }) {
        my @wells = @{ $relation->{wells} };

        my $ph_name = $relation->{miseq_dets}->{name} . '_FP';
        print "\nLooking for placeholder - $ph_name";
        my $placeholder_fp = $lims2_model->schema->resultset('Plate')->find({ name => $ph_name });

        print "\n~~~~~~~~~~~~~~~~~~~~~~~~~ Linking $piq_key -> $relation->{miseq_dets}->{name} ~~~~~~~~~~~~~~~~~~~~\n";
        foreach my $process (@wells) {
            if ($placeholder_fp) {
                print "\nChecking for Placeholder well";
                my $placeholder_well = $lims2_model->schema->resultset('Well')->find({
                    name => $process->{output_wells}[0]->{well_name},
                    plate_id => $placeholder_fp->id
                });
                if ($placeholder_well) {
                    my @ph_processes = map { $_->id } $placeholder_well->child_processes;
                    foreach my $ph_process (@ph_processes) {
                        print "\nDeleting place-holder process: $ph_process";
                        $lims2_model->delete_process({ id => $ph_process });
                    }
                    print "\nDeleting place-holder well: $placeholder_fp $process->{output_wells}[0]->{well_name}";
                    $lims2_model->delete_well({ id => $placeholder_well->id });
                }
            }
            my $process_rs = $lims2_model->create_process($process);
        }
    }
}
print "\nUpdating 2nd QC Miseq Experiments\n";
my $seen;
my @uniq_piqs = grep { !$seen->{$_->{name}} ++ } @eps_piqs;

foreach my $piq_qc (@uniq_piqs) {
    my $piq_exp_rs = $lims2_model->schema->resultset('MiseqExperiment')->search({
        name => { -like => $piq_qc->{name} . '%' },
    });
    my $piq_rs = $lims2_model->schema->resultset('Plate')->find({ name => $piq_qc->{piq} });
    while (my $piq_exp = $piq_exp_rs->next) {
        $lims2_model->update_miseq_experiment({ 
            id              => $piq_exp->id,
            experiment_id   => $piq_qc->{experiment_id},
            parent_plate_id => $piq_rs->id,
            gene            => $piq_qc->{gene},
        });
    }
}

if ($dump) {
    print "\n~~~~~~~~~~~~~~~~~ PIQ parental info ~~~~~~~~~~~~~~~~~\n";

    print Dumper $piq_data;

    print "\n~~~~~~~~~~~~~~~~~ Miseq child info  ~~~~~~~~~~~~~~~~~\n";

    print Dumper $miseq_data;
} else { 
    open(my $pfh, '>', 'piq_data.txt');
    print $pfh Dumper $piq_data;
    close $pfh;

    open(my $mfh, '>', 'miseq_data.txt');
    print $mfh Dumper $miseq_data;
    close $mfh;
}

print "\n~~~~~~~~~~~~~~~~~ Unknown Entities  ~~~~~~~~~~~~~~~~~\n";

if ($unknown_entities) {
    print Dumper $unknown_entities;
} else {
    print "\n\t\t\tNil\n\n";
}

1;
