#!/usr/bin/env perl
use strict;
use warnings;
use File::Slurp qw( read_dir );
use File::Spec::Functions qw( catfile );
use Path::Class;
use Data::Dumper;
use LIMS2::Model;
use Try::Tiny;

my $warehouse = dir("/warehouse/team229_wh01/lims2_managed_sequencing_data");
my @dirs = grep { -d } map { catfile $warehouse, $_ } read_dir $warehouse;
my @backups;
foreach my $project (@dirs) {
    my @split_sub = split('/',$project);
    my $seq->{name} = $split_sub[-1];
    print $project . "\n";
    my $seq_dir = dir("/warehouse/team229_wh01/lims2_managed_sequencing_data", $seq->{name});
    my @sub_dirs = grep { -d } map { catfile $seq_dir, $_ } read_dir $seq_dir;
    foreach my $subs(@sub_dirs) {
        print "Found backup: " . $subs . "\n";
        my @split_sub = split('/',$subs);
        my $epoch_timestamp = (stat($subs))[9];
        my @date=localtime($epoch_timestamp);
        $date[5] +=1900;
        $date[4] +=1;
        $seq->{$split_sub[-1]} = "$date[5]-$date[4]-$date[3] $date[2]:$date[1]:$date[0]";
    }
    if (keys %{$seq} > 1) {
        push(@backups, $seq);
    }
}
print "\n-------------------------------------------\n\n";
print Dumper @backups;
print "\n-------------------------------------------\n\n";

my $model = LIMS2::Model->new( user => 'lims2' );
foreach my $record (@backups) {
    my @attrs = keys %{$record};
    foreach my $attr (@attrs) {
        if ($attr ne "name") {
            try{
                my $seq_rs = $model->schema->resultset('SequencingProject')->find({
                    name => $record->{name},
                })->{_column_data};
                print "(" . $seq_rs->{id} . ", \'" . $attr . "\', \'" . $record->{$attr} . "\'),\n";
            } catch {
                #print $record->{name} . " not found.\n";
            };
        }
    }
} 
