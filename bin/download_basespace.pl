#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use LIMS2::Model::Util::BaseSpace;

my @samples = ();
my $path = q/./;
GetOptions(
    'sample=i@' => \@samples,
    'path=s'    => \$path)
    or die 'Error in command line arguments';

my $api  = LIMS2::Model::Util::BaseSpace->new;
foreach my $sample_id ( @samples ) {
my $sample = $api->sample($sample_id);
    foreach my $file ( $sample->files ) {
        printf "File id: %s\nName: %s\nE-tag: %s\n",
            $file->id, $file->name, $file->etag;
       $file->download($path);
    }
}

