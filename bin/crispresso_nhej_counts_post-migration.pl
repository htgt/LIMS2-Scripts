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
use LIMS2::Model::Util::ImportCrispressoQC qw( get_crispr );
use List::Compare::Functional qw( get_intersection );
use Scalar::Util qw(looks_like_number);

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
    if ($common_read[2] eq 'True' ) {
        $fs_check = ($common_read[5] + $common_read[6]) % 3;
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

sub header_hash {
    my ( $header, @expected_titles ) = @_;
    my @titles = split( /\t/, lc $header );
    my @intersection = get_intersection( [ \@titles, \@expected_titles ] );
    my %head;
    %head = map { lc $titles[$_] => $_ }
        0 .. $#titles;    #creates a hash that has all the elements of the array as keys and their index as values

#check if the length of the intersection of the full array of titles is equal to the length of the array of expected titles
#This checks that all the requested elements were found within the header of the file
    if (scalar(@intersection) == scalar(@expected_titles)) {
        return %head;
    }
    else {
        warn "Miseq Alleles Frequency file does not hold all the required headers";
        return;
    }
    return 1;
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

    foreach my $exp (@exps) {
        my $quant = find_file($base, $i, $exp, "Quantification_of_editing_frequency.txt");
        if ($quant) {
            my $fh;
            open ($fh, '<:encoding(UTF-8)', $quant) or die "$!";
            chomp(my @lines  = <$fh>);
            close $fh;

            my %params;
            my %commands = (
                'NHEJ'           => 'nhej_reads',
                'HDR'            => 'hdr_reads',
                'Mixed HDR-NHEJ' => 'mixed_reads',
                'Total Aligned'  => 'total_reads',
            );
            while ( @lines) {
                my $line =  shift @lines;
                my ( $type, $number ) = $line =~ m/^[\s\-]* #ignore whitespace, dashes at start of line
                    ([\w\-\s]+) #then grab the type e.g. "HDR", "Mixed HDR-NHEJ", etc
                    :(\d+) #finally the number of reads
                    /xms;
                next unless $type && $number;
                if ( exists $commands{$type} ) {
                    $params{ $commands{$type} } = $number;
                }

            }


            $experiments->{$exp}->{$i} = {
                        nhej_reads      => $params{nhej_reads},
                        total_reads     => $params{total_reads},
                        hdr_reads       => $params{hdr_reads},
                        mixed_reads     => $params{mixed_reads},
                    };
        }
        my $read = find_file($base, $i, $exp, "Alleles_frequency_table.txt");

        if ($read) {
            my $fh;
            open ($fh, '<:encoding(UTF-8)', $read) or die "$!";
            chomp(my @lines  = <$fh>);
            close $fh;
            $experiments->{$exp}->{$i}->{classification} = 'Not Called';
            $experiments->{$exp}->{$i}->{frameshifted} = 0;

            if (scalar @lines > 3) {
                my @mixed_read = split(/\t/, $lines[3]);
                my $mixed_check = $mixed_read[$#mixed_read];
                if ($mixed_check >= 5) {
                    $experiments->{$exp}->{$i}->{frameshifted} = 0;
                    $experiments->{$exp}->{$i}->{classification} = 'Mixed';
                } else {
                    my @first_most_common = split(/\t/, $lines[1]);
                    my @second_most_common = split(/\t/, $lines[2]);

                    my $fs_check = frameshift_check($experiments, @first_most_common) + frameshift_check($experiments, @second_most_common);
                    if ($fs_check != 0) {
                        $experiments->{$exp}->{$i}->{classification} = 'Not Called';
                        $experiments->{$exp}->{$i}->{frameshifted} = 1;
                    }
                }
            }


            my $limit = 10;
            my $counter = 0;
            my %histogram;
            my $header = shift(@lines);    #grab the header line that holds the titles of the columns
            my @expected_titles = ( 'aligned_sequence', 'nhej', 'unmodified', 'hdr', 'n_deleted', 'n_inserted', 'n_mutated', '#reads' );
            my %head = header_hash( $header, @expected_titles );

            while (@lines) {
                my $line = shift @lines;
                my @elements = split /\t/ , $line;
                my $indel;
                if (looks_like_number($elements[$head{n_inserted}]) and looks_like_number($elements[$head{n_deleted}])) {
                    $indel = int ($elements[$head{n_inserted}]) - int ($elements[$head{n_deleted}]);
                }
                else {
                    print "Headers are not working as intended";
                }
                if ($indel) {
                    $histogram{$indel}+=$elements[$head{'#reads'}];
                } 
                if ( $counter < $limit ) {
                    $counter++;
                    my $row = {
                        aligned_sequence         => $elements[ $head{aligned_sequence} ],
                        nhej                     => lc $elements[ $head{nhej} ],
                        unmodified               => lc $elements[ $head{unmodified} ],
                        hdr                      => lc $elements[ $head{hdr} ],
                        n_deleted                => int $elements[ $head{n_deleted} ],
                        n_inserted               => int $elements[ $head{n_inserted} ],
                        n_mutated                => int $elements[ $head{n_mutated} ],
                        n_reads                  => int $elements[ $head{'#reads'} ],
                    };
                    if ($head{phred_quality}){
                        $row->{quality_score} = $elements[ $head{phred_quality}];

                    }
                    if ($head{reference_sequence}){
                        $row->{reference_sequence} = $elements[ $head{reference_sequence}];

                    }

                    push ( @{$experiments->{$exp}->{$i}->{allele_frequencies}}, $row );
                }
            }
            $experiments->{$exp}->{$i}->{histogram} = \%histogram;
        }
        my $jobout = find_file($base, $i, $exp, "job.out");
        if ($jobout) {
            my $job = get_crispr($jobout);
            $experiments->{$exp}->{$i}->{jobout} = $job;
        }
    }
}

my $result;
foreach my $exp (keys %{$experiments}) {
    my $nhej = 0;
    my $total = 0;
    foreach my $index (keys %{$experiments->{$exp}}) {
        my $nhej_well = $experiments->{$exp}->{$index}->{nhej_reads} || 0;
        my $total_well = $experiments->{$exp}->{$index}->{total_reads} || 0;
        $nhej += $nhej_well;
        $total += $total_well;
    }
    unless ($total) {
        next;
    }
    my $target = sprintf("%.2f", ($nhej / $total) * 100);
    $result->{$exp} = {
        nhej    => $nhej,
        total   => $total,
        efficiency => $target,
    };
    print "Experiment: " . $exp . ", NHEJ: " . $nhej . ", Total: " . $total . ", Eff: " . $target . "%\n";
}

if ($summary) { #One time use code
    my $csv = Text::CSV->new({binary => 1, eol => "\n"}) or die "Cannot use CSV: ".Text::CSV->error_diag ();

    my $old_file = $base . 'summary.csv';
    my $new_file = $old_file . '.tmp';

    open my $in, "<:encoding(utf8)", $old_file or die "$old_file: $!";
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
                        nhej_reads      => $result->{$exp}->{nhej} || '0',
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
        my @wells = keys %{$experiments->{$exp}};
        for (my $well = 1; $well < 385; $well++) {
            my $well_rs = $model->schema->resultset('Well')->find({ plate_id => $plate_rs->{id}, name => $well_names[$well - 1] });
            if ($well_rs) {
                $well_rs = $well_rs->as_hash;
            }
            else {
                print "Plate: " . $plate_rs->{id} . ", Well: " . $well_names[$well - 1] . " - Failed to retrieve well \n";
            next;
            }
            my $well_exp = $model->schema->resultset('MiseqWellExperiment')->find({ well_id => $well_rs->{id}, miseq_exp_id => $exp_check->{id} });
            unless ( grep( /^$well$/, @wells ) ) {
                if ($well_exp){
                    $well_exp = $well_exp->as_hash;
                    $model->schema->resultset('MiseqAllelesFrequency')->search( { miseq_well_experiment_id => $well_exp->{id} } )->delete_all;
                    $model->schema->resultset('IndelHistogram')->search({ miseq_well_experiment_id => $well_exp->{id}} )->delete_all;
                    $model->schema->resultset('CrispressoSubmission')->search({ id => $well_exp->{id}})->delete_all;
                    $model->schema->resultset('MiseqWellExperiment')->search( { id => $well_exp->{id} } )->delete_all;
                    print "Deleted Miseq Well Exp ID: " . $well_exp->{id} . "\n";
                }
            }
            else {
                if ($experiments->{$exp}->{$well}->{total_reads}) {
                    print "Attempt Well: " . $well . "\n";
                    unless ($well_exp) {
                        $model->schema->txn_do( sub {
                            try {
                                $well_exp = $model->create_miseq_well_experiment({
                                    well_id         => $well_rs->{id},
                                    miseq_exp_id    => $exp_check->{id},
                                    classification  => $experiments->{$exp}->{$well}->{classification},
                                    frameshifted    => $experiments->{$exp}->{$well}->{frameshifted},
                                    total_reads     => $experiments->{$exp}->{$well}->{total_reads} || 0, #UNTESTED
                                    nhej_reads      => $experiments->{$exp}->{$well}->{nhej_reads} || 0, #UNTESTED
                                    hdr_reads       => $experiments->{$exp}->{$well}->{hdr_reads} || 0, #UNTESTED
                                    mixed_reads     => $experiments->{$exp}->{$well}->{mixed_reads} || 0, #UNTESTED
                                })->as_hash;
                            print "Created Miseq Well Exp ID: " . $well_exp->{id} . "\n"
                            }
                            catch {
                                warn "Could not create well record for " . $well_rs->{id} . ": $_";
                            };
                        });
                    }
                    else {
                    $well_exp = $well_exp->as_hash;
                    $model->schema->resultset('MiseqAllelesFrequency')->search( { miseq_well_experiment_id => $well_exp->{id} } )->delete_all;
                    $model->schema->resultset('IndelHistogram')->search({ miseq_well_experiment_id => $well_exp->{id}} )->delete_all;
                    $model->schema->resultset('CrispressoSubmission')->search({ id => $well_exp->{id}})->delete_all;

                    $model->schema->txn_do( sub {
                        try {
                            $model->update_miseq_well_experiment({
                                id              => $well_exp->{id},
                                classification  => $well_exp->{classification} || $experiments->{$exp}->{$well}->{classification},
                                frameshifted    => $experiments->{$exp}->{$well}->{frameshifted},
                                nhej_reads      => $experiments->{$exp}->{$well}->{nhej_reads} || 0, #UNTESTED
                                total_reads     => $experiments->{$exp}->{$well}->{total_reads} || 0, #UNTESTED
                                hdr_reads       => $experiments->{$exp}->{$well}->{hdr_reads} || 0, #UNTESTED
                                mixed_reads     => $experiments->{$exp}->{$well}->{mixed_reads} || 0, #UNTESTED
                            });
                        print "Updated Miseq Well Exp ID: " . $well_exp->{id} . " Frameshifted:" . $experiments->{$exp}->{$well}->{frameshifted} . " Well: " . $well . "\n";
                        }
                        catch {
                            warn "Could not update well record for " . $well_exp->{id} . ": $_";
                        };
                    });
                }
                my @alleles = @{$experiments->{$exp}->{$well}->{allele_frequencies}};
                if (@alleles) {
                    foreach my $freq (@alleles){
                        $freq->{miseq_well_experiment_id} = $well_exp->{id};
                        $model->schema->txn_do(
                            sub {
                                try {
                                    $model->create_miseq_alleles_frequency($freq);
                                }
                                catch {
                                    warn "Error creating entry";
                                    $model->schema->txn_rollback;
                                };
                            }
                        );
                    }
                }

                my $histo = $experiments->{$exp}->{$well}->{histogram};
                if ($histo) {
                    foreach my $key (keys %{$histo}) {
                        my $row = {
                            miseq_well_experiment_id    =>  $well_exp->{id},
                            indel_size                  =>  $key,
                            frequency                   =>  $histo->{$key},
                        };
                        $model->schema->txn_do(
                            sub {
                                try {
                                    $model->create_indel_histogram($row);
                                }
                                catch {
                                    warn "Error creating indel histogram entry";
                                    $model->schema->txn_rollback;
                                };
                            }
                        );
                    }
                }

                my $jobout = $experiments->{$exp}->{$well}->{jobout};
                if ($jobout) {
                    $jobout->{id} = $well_exp->{id};
                    $model->schema->txn_do(
                        sub {
                            try {
                                $model->create_crispr_submission($jobout);
                            }
                            catch {
                                warn "Error creating cripr submission entry";
                                $model->schema->txn_rollback;
                            };
                        }
                    );
                }
            }
            else {
                if ($well_exp->{id}){
                    $model->schema->resultset('MiseqAllelesFrequency')->search( { miseq_well_experiment_id => $well_exp->{id} } )->delete_all;
                    $model->schema->resultset('IndelHistogram')->search({ miseq_well_experiment_id => $well_exp->{id}} )->delete_all;
                    $model->schema->resultset('CrispressoSubmission')->search({ id => $well_exp->{id}})->delete_all;
                    $model->schema->resultset('MiseqWellExperiment')->search( { id => $well_exp->{id} } )->delete_all;
                    print "Deleted Miseq Well Exp ID: " . $well_exp->{id} . "\n";
                }
            }
        }
    }
    }
    print "Finished.";
}
