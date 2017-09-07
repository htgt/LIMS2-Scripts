#! /usr/bin/perl

use strict;
use warnings;

use feature qw(say);
use Getopt::Long;
use Data::Dumper;
use File::Find;
use File::Copy qw(move);
use Text::CSV;
use LIMS2::Model;
use Try::Tiny;

GetOptions(
    'project=s' => \my $project,
    'db_update' => \my $db_update,
    'summary'   => \my $summary,
);

sub find_children {
    my ( $base, $reg ) = @_;
    my $fh;
    opendir ($fh, $base);
    my @files = grep {/$reg/} readdir $fh;
    closedir $fh;
    return @files;
}

sub find_file {
    my ( $base, $index, $exp, $name ) = @_;

    my $nhej_files = [];
    my $wanted = sub { _wanted( $nhej_files, $name ) };
    my $dir = $base . "S" . $index . "_exp" . $exp;
    find($wanted, $dir);

    return @$nhej_files[0];
}

sub _wanted {
    return if ! -e;
    my ( $nhej_files, $file_name ) = @_;

    push(@$nhej_files, $File::Find::name) if $File::Find::name=~ /$file_name/;

    return;
}

sub read_file_lines {
    my ( $fh ) = @_;

    my @data;
    while (my $row = <$fh>) {
        chomp $row;
        push(@data, join(',', split(/\t/,$row)));
    }

    return @data;
}

sub read_columns {
    my ( $model, $csv, $fh ) = @_;

    my $overview;

    while ( my $row = $csv->getline($fh)) {
        next if $. < 2;
        my @genes;
        push @genes, $row->[1];
        $overview->{$row->[0]} = \@genes;
    }

    return $overview;
}

sub file_handling {
    my $file_name = shift;

    my $fh;
    open ($fh, '<:encoding(UTF-8)', $file_name) or die "$!";
    my @lines = read_file_lines($fh);
    close $fh;
    
    return \@lines;
}

sub frameshift_check {
    my ($experiments, @common_read) = @_; 
    my $fs_check = 0;
    if ($common_read[1] eq 'True' ) {
        $fs_check = ($common_read[4] + $common_read[5]) % 3;
    }
    return $fs_check;
}

sub well_builder {
    my ($mod, @well_names) = @_;
    
    foreach my $number (1..12) {
        my $well_num = $number + $mod->{mod};
        foreach my $letter ( @{$mod->{letters}} ) {
            my $well = sprintf("%s%02d", $letter, $well_num);
            push (@well_names, $well);
        }
    }

    return @well_names;
}

my $rna_seq = $ENV{LIMS2_RNA_SEQ} || "/warehouse/team229_wh01/lims2_managed_miseq_data/";
my $base = $rna_seq . $project . '/';
my $experiments;

for (my $i = 1; $i < 385; $i++) {
    my $reg = "S" . $i . "_exp[A-Za-z0-9_]+";
    my @files = find_children($base, $reg);
    my @exps;
    foreach my $file (@files) {
        my @matches = ($file =~ /S\d+_exp([A-Za-z0-9_]+)/g);
        foreach my $match (@matches) {
            push (@exps,$match);
        }
    }

    @exps = sort @exps;
    my @selection;
    my $percentages;
    my $classes;

    foreach my $exp (@exps) {
        my $quant = find_file($base, $i, $exp, "Quantification_of_editing_frequency.txt");

        if ($quant) {
            my @lines = @{file_handling($quant)};
            my @nhej = ($lines[2] =~ qr/^,- NHEJ:(\d+)/);
            my @total = ($lines[6] =~ qr/^Total Aligned:(\d+)/);

            $experiments->{$exp}->{sprintf("%02d", $i)} = {
                nhej    => $nhej[0],
                total   => $total[0],
            };
        }

        my $read = find_file($base, $i, $exp, "Alleles_frequency_table.txt");

        if ($read) {
            my @lines = @{file_handling($read)};
            if (scalar @lines > 3) {
                my @mixed_read = split(/,/, $lines[3]);
                my $mixed_check = $mixed_read[$#mixed_read];
                if ($mixed_check >= 5) {
                    $experiments->{$exp}->{sprintf("%02d", $i)}->{frameshifted} = 0;
                    $experiments->{$exp}->{sprintf("%02d", $i)}->{classification} = 'Mixed';
                } else {
                    my @first_most_common = split(/,/, $lines[1]); 
                    my @second_most_common = split(/,/, $lines[2]);
                   
                    my $fs_check = frameshift_check($experiments, @first_most_common) + frameshift_check($experiments, @second_most_common);
                    
                    if ($fs_check != 0) {
                        $experiments->{$exp}->{sprintf("%02d", $i)}->{classification} = 'Not Called';
                        $experiments->{$exp}->{sprintf("%02d", $i)}->{frameshifted} = 1;
                    }
                }
            }
        }
    }
}

my $result;

foreach my $exp (keys %{$experiments}) {
    my $nhej = 0;
    my $total = 0;
    foreach my $index (keys %{$experiments->{$exp}}) {
        $nhej += $experiments->{$exp}->{$index}->{nhej};
        $total += $experiments->{$exp}->{$index}->{total};
    }
    my $target = sprintf("%.2f", ($nhej / $total) * 100);
    $result->{$exp} = {
        nhej    => $nhej,
        total   => $total,
        efficiency => $target,
    };
}

if ($summary) { #One time use code
    my $csv = Text::CSV->new({binary => 1, eol => "\n"}) or die "Cannot use CSV: ".Text::CSV->error_diag ();

    my $old_file = $base . 'summary.csv';
    my $new_file = $old_file . '.tmp';

    open my $in, $old_file or die "$old_file: $!";
    open my $out, '>', $new_file or die "$new_file: $!";
    
    my $header = $csv->getline($in);

    if (scalar @$header != 9) {
        splice @$header, 6, 0, "NHEJ";
        splice @$header, 7, 0, "Total";
        $csv->print($out, $header);
        while (my $row = $csv->getline($in)) {
            splice @$row, 6, 0, $result->{@$row[0]}->{nhej};
            splice @$row, 7, 0, $result->{@$row[0]}->{total};
            $csv->print($out, $row);
        }
    } else {
        $csv->print($out, $header);
        while (my $row = $csv->getline($in)) {
            $csv->print($out, $row);
        }
    }

    close $in;
    close $out;

    move($new_file, $old_file) or die "Can't move: $!";
}

if ($db_update) {
    #Update miseq_exp table
    my $model = LIMS2::Model->new({ user => 'tasks' });

    my $csv = Text::CSV->new({ binary => 1 }) or die "Can't use CSV: " . Text::CSV->error_diag();
    open my $fh, '<:encoding(UTF-8)', $base . '/summary.csv' or die "Can't open CSV: $!";
    my $ov = read_columns($model, $csv, $fh);
    close $fh;

    my $plate_rs = $model->schema->resultset('Plate')->find({ name => $project })->as_hash;
    my $proj_rs = $model->schema->resultset('MiseqPlate')->find({ plate_id => $plate_rs->{id} })->as_hash;

    my @well_names;
    my $quads = {
        '0' => {
            mod     => 0,
            letters => ['A','B','C','D','E','F','G','H'],
        },
        '1' => {
            mod     => 12,
            letters => ['A','B','C','D','E','F','G','H'],
        },
        '2' => {
            mod     => 0,
            letters => ['I','J','K','L','M','N','O','P'], 
        },
        '3' => {
            mod     => 12,
            letters => ['I','J','K','L','M','N','O','P'], 
        }
    };

    for (my $ind = 0; $ind < 4; $ind++) {
        @well_names = well_builder($quads->{$ind}, @well_names);
    }

    foreach my $exp (keys %{$result}) {
        my $exp_check = $model->schema->resultset('MiseqExperiment')->find({ miseq_id => $proj_rs->{id}, name => $exp });
        unless ($exp_check) {
            $model->schema->txn_do( sub {
                try {
                    $model->create_miseq_experiment({
                        miseq_id        => $proj_rs->{id},
                        name            => $exp,
                        gene            => (split(/_/,$ov->{$exp}[0]))[0],
                        mutation_reads  => $result->{$exp}->{nhej} || '0',
                        total_reads     => $result->{$exp}->{total} || '1',
                    });
                    print "Inserted Miseq ID: " . $proj_rs->{id} . " Experiment: " . $exp . "\n";
                }
                catch {
                    warn "Could not create record for " . $proj_rs->{id} . ": $_";
                    $model->schema->txn_rollback;
                };
            });
            $exp_check = $model->schema->resultset('MiseqExperiment')->find({ miseq_id => $proj_rs->{id}, name => $exp });
        }

        $exp_check = $exp_check->as_hash;
        foreach my $well (keys %{$experiments->{$exp}}) {
            if (defined $experiments->{$exp}->{$well}->{frameshifted}) {
                print "Attempt Well: " . $well . "\n";
                my $well_rs = $model->schema->resultset('Well')->find({ plate_id => $plate_rs->{id}, name => $well_names[$well - 1] })->as_hash;
                my $well_exp = $model->schema->resultset('MiseqWellExperiment')->find({ well_id => $well_rs->{id}, miseq_exp_id => $exp_check->{id} });

                if ($well_exp) {
                    $well_exp = $well_exp->as_hash;

                    $model->schema->txn_do( sub {
                        try {
                            $model->update_miseq_well_experiment({
                                id              => $well_exp->{id},
                                classification  => $well_exp->{classification} || $experiments->{$exp}->{$well}->{classification},
                                frameshifted    => $experiments->{$exp}->{$well}->{frameshifted},
                            });
                            print "Updated Miseq Well Exp ID: " . $well_exp->{id} . " Frameshifted:" . $experiments->{$exp}->{$well}->{frameshifted} . " Well: " . $well . "\n";
                        }
                        catch {
                            warn "Could not update well record for " . $well_exp->{id} . ": $_";
                        };
                    });
                } else {
                    $model->schema->txn_do( sub {
                        try {
                            $model->create_miseq_well_experiment({
                                well_id         => $well_rs->{id},
                                miseq_exp_id    => $exp_check->{id},
                                classification  => $experiments->{$exp}->{$well}->{classification},
                                frameshifted    => $experiments->{$exp}->{$well}->{frameshifted},
                            });
                            $well_exp = $model->schema->resultset('MiseqWellExperiment')->find({ well_id => $well_rs->{id}, miseq_exp_id => $exp_check->{id} })->as_hash;
                            print "Created Miseq Well Exp ID: " . $well_exp->{id} . " Frameshifted:" . $experiments->{$exp}->{$well}->{frameshifted} . "\n";
                        }
                        catch {
                            warn "Could not create well record for " . $well_rs->{id} . ": $_";
                        };
                    });
                }
            }
        }
    }
}
