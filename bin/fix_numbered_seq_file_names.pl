#! /usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;
use Path::Class;

my ($dir_name,$project,$primer) = @ARGV;

# Build array of well names in order a01,a02,a03....b01,b02 and so on
my @well_names;
foreach my $letter ( qw(a b c d e f g h)){
	foreach my $number (1..12){
		my $well = sprintf("%s%02d",$letter,$number);
		push @well_names, $well;
	}
}

my $dir = dir($dir_name);

my @files = $dir->children;
foreach my $data_file (@files){
    my $old_name = $data_file->basename;

    # Files named by sample number instead of well name,
    # e.g. HFP0015_A_P19F#10.seq
    my ($number) = ( $old_name =~ /\#(\d*)/g );

    # When all the samples have the same name the files
    # for the first sample are named with no number e.g. HFP0015_A_P19F.seq
    if(! defined $number){
        my $a01_file_name = $project."_".$primer;
        if ($old_name =~ /$a01_file_name\./){
            $number = 1;
        }
    }

    my ($suffix) = ( $old_name =~ /\.([^\.]*)$/ );
    unless($number and $suffix){
    	say "Error: cannot indentify well number and suffix in file name $old_name";
    	next;
    }
    my $well = $well_names[$number-1];
    my $new_name = $project.$well.".p1k".$primer.".$suffix";
    my $new_path = $dir->file($new_name);
    say "moving $data_file to $new_path";
    $data_file->move_to($new_path);
}