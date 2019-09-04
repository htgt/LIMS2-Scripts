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

my @files = grep (/.*\.json$/, @paths);
my $json = JSON->new;

my $name_regex = '^\ ?(\d+)\)?\ ?([A-z0-9]+)\ *([A-z0-9]+)\ ?\ *([A-Z][\dOI]{2})\ +(HET|HOM|WT).*$';
my $tt_regex = '(Plate|PLATE|plate)\_?\ ?(\d+)';


# 2)HUPFP0033_2 RNF43_3 D06 WT  || 1)HUPFP0033_2 RNF43_3 D05 HET BN
#(2)(HUPFP0033_2)      (D06)    ||(1)(HUPFP0033_2)      (D05)

my @wells = generate_piq_well_names(1);

my $unknown_entities;
my $piq_data;
my $miseq_data;

foreach my $json_filepath (@files) {
    my $miseq_check = 1;
    my $full_path = $dir . $json_filepath;
    my $json_text = do {
        open(my $json_fh, "<:encoding(UTF-8)", $full_path)
            or die("Can't open $json_filepath: $!\n");
        local $/;
        my $json_data = <$json_fh>;
        close $json_fh;
        $json_data
    };

    my $card_data = $json->decode($json_text);

    my $piq_name = $card_data->{name};
    my $card_desc = $card_data->{desc};
    my @miseq_captures = $card_desc =~ /(Miseq.\d+)/gmi;
    my $miseq_name;
    if (@miseq_captures) {
        $miseq_name = $miseq_captures[0];
    }
    $piq_data->{$piq_name} = {
        name            => $piq_name,
        created_by      => 'pk8@sanger.ac.uk',
        species         => 'Human',
        created_at      => strftime("%Y-%m-%dT%H:%M:%S", localtime(time)),
        type            => 'PIQ',
        miseq           => $miseq_name,
    };
    
    
    my $arr_num = 0;
    print Dumper $piq_name;

    my @parent_data = @{ $card_data->{checklists}[$arr_num]->{checkItems} };
    my @piq_wells;
    my @miseq_processes;
    foreach my $parent (@parent_data) {
        my $name_str = $parent->{name};
        my ($piq_inc, $plate_name, $gene, $well_name, $primary_call) = $name_str =~ m/$name_regex/;
        if ($plate_name) {
            $plate_name = remove_well_leading_zero($plate_name);
        } else {
            push ( @{$unknown_entities->{$piq_name}->{names}}, $name_str );
            $unknown_entities->{$piq_name}->{reason} = 'Did not return a plate name';
            next;
        }
        my $plate_rs = $lims2_model->schema->resultset('Plate')->search({ name => $plate_name });
        if ($plate_rs->count == 0) {
            $plate_rs = $lims2_model->schema->resultset('Plate')->search({ name => $gene });
        }
        my $tt_plate;
        if ($plate_rs->count == 0) {
            my ($plate_tag, $tt_id) = $name_str =~ m/$tt_regex/;
            if ($tt_id){
                $tt_plate = 'HUPFPPlate' . $tt_id;
                $plate_rs = $lims2_model->schema->resultset('Plate')->search({ name => $tt_plate });
            }   
        }
        if ($plate_rs->count == 0) {
            if ($tt_plate) {
                $unknown_entities->{$piq_name}->{$plate_name} = { relation => $name_str, reason => "Plate $plate_name & $tt_plate not found." };
            } else {
                $unknown_entities->{$piq_name}->{$plate_name} = { relation => $name_str, reason => "Plate $plate_name not found." };
            }
            next;
        }
        $well_name =~ s/O/0/g; 
        $well_name =~ s/I/1/g; 
        my $well_rs = $lims2_model->schema->resultset('Well')->search({ 
            name        => $well_name,
            plate_id    => { -in => $plate_rs->get_column('id')->as_query },
        });
        if ($well_rs->count == 0) {
            $unknown_entities->{$piq_name}->{$plate_name} = { well => $well_name, relation => $name_str, reason => 'Well not found' };
            next; 
        }

        my $gene_symbol = (split '_', $gene)[0];
        my @miseq_exps = map { { miseq => $_->miseq_plate->{plate}, exp_name => $_->name, exp_id => $_->id } } $lims2_model->schema->resultset('MiseqExperiment')->search({
            name => {-like => $piq_name . '%' },
            gene => {-like => $gene_symbol . '%'},
        });
        my $piq_well_name = $wells[$piq_inc - 1];
        my $well_data = {
            well_name       => $piq_well_name,
            process_type    => 'dist_qc',
            parent_well     => $well_rs->first->name,
            parent_plate    => $plate_rs->first->name,
            primary_call    => $primary_call,
            gene            => $gene,
            secondary_exp   => @miseq_exps,
        };
        push (@piq_wells, $well_data);
    }

    @piq_wells = sort { $a->{well_name} cmp $b->{well_name} } @piq_wells;
    $piq_data->{$piq_name}->{wells} = \@piq_wells;
    push @{$miseq_data->{$piq_name}} , {
        piq         => $piq_name,
        wells       => \@miseq_processes,
    };

}


if ($dump) {
    # print "\n~~~~~~~~~~~~~~~~~ PIQ parental info ~~~~~~~~~~~~~~~~~\n";

    #print Dumper $piq_data;

    #print "\n~~~~~~~~~~~~~~~~~ Miseq child info  ~~~~~~~~~~~~~~~~~\n";

    #print Dumper $miseq_data;
} else { 
    open(my $pfh, '>', 'piq_data.txt');
    print $pfh Dumper $piq_data;
    close $pfh;

    open(my $mfh, '>', 'miseq_data.txt');
    print $mfh Dumper $miseq_data;
    close $mfh;
}

#print "\n~~~~~~~~~~~~~~~~~ Unknown Entities  ~~~~~~~~~~~~~~~~~\n";

if ($unknown_entities) {
    # print Dumper $unknown_entities;
} else {
    print "\n\t\t\tNil\n\n";
}

my $csv = Text::CSV->new({binary => 1, eol => $/ })
    or die "Failed to create a CSV handle: $!";

open(my $mfh, '>', 'unknown_data.txt');
print $mfh Dumper $unknown_entities;
close $mfh;

my @uheaders = qw(PIQ_plate FP_plate well name_string reason);
open my $ufh, ">:encoding(utf8)", 'unknown_entries.csv' or die "failed to create unknown.csv: $!";
$csv->print($ufh, \@uheaders);
foreach my $piq (keys %{$unknown_entities}) {
    my $piqh = $unknown_entities->{$piq};
    foreach my $str (@{$piqh->{names}}) {
        my @row = ($piq,'','',$str,$piqh->{reason});
        $csv->print($ufh, \@row);
    }
    delete $piqh->{names};
    delete $piqh->{reason};
    foreach my $plate (keys %{$piqh}) {
        my @row = ($piq, $plate,'',$piqh->{$plate}->{relation}, $piqh->{$plate}->{reason});
        if ($piqh->{$plate}->{well}) {
            @row = ($piq,$plate,$piqh->{$plate}->{well},$piqh->{$plate}->{relation}, $piqh->{$plate}->{reason});
        }
        $csv->print($ufh, \@row);
    }
}
close $ufh;

my $filename = 'secondary_qc_plates.csv';
my @headers = qw(PIQ_Plate PIQ_Well Gene Miseq_Plate Parent_Plate Parent_Well Primary_Call Miseq_Experiment Miseq_Exp_Plate Miseq_Exp_ID );

open my $fh, ">:encoding(utf8)", $filename or die "failed to create $filename: $!";
my @keys = keys %{ $piq_data };
@keys = sort @keys;
$csv->print($fh, \@headers);
foreach my $piq (@keys) {
    foreach my $well (@{ $piq_data->{$piq}->{wells} }) {
        my @row = ($piq, $well->{well_name}, $well->{gene}, $piq_data->{$piq}->{miseq}, $well->{parent_plate}, $well->{parent_well}, $well->{primary_call});
        my @exps = $well->{secondary_exp};
        if (@exps) {
            foreach my $exp (@exps) {
                my @exp_row = @row;
                push @exp_row, $exp->{exp_name}, $exp->{miseq}, $exp->{exp_id};
                $csv->print($fh, \@exp_row);
            }
        } else {
            $csv->print($fh, \@row);
        }
    }
}
close $fh;

1;
