#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use LIMS2::Model::Util::BaseSpace;

my $walkup;
my $path = q/./;
GetOptions(
    'walkup_name=s' => \$walkup,
    'path=s'    => \$path)
    or die 'Error in command line arguments';

my $api  = LIMS2::Model::Util::BaseSpace->new;

my $undetermined_samples = $api->get("projects/34207198/samples");
my $walkup_ids;
foreach my $project (@{ $undetermined_samples->{Items} }) {
    $walkup_ids->{ $project->{ExperimentName} } = $project->{Id};
}

my $undetermined_sample = $api->sample($walkup_ids->{$walkup});
foreach my $file ( $undetermined_sample->files ) {
    printf "File id: %s\nName: %s\nE-tag: %s\n",
        $file->id, $file->name, $file->etag;
    $file->download($path);
}

1;
