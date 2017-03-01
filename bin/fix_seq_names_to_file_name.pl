#! /usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;
use Path::Class;

my @dir_names = @ARGV;

foreach my $dir_name (@dir_names){
	my $dir = dir($dir_name);

	while(my $data_file = $dir->next){
	    my $old_name = $data_file->basename;
	    my ($base, $suffix) = ( $old_name =~ /(.*)\.([^\.]*)$/ );
	    if( $suffix eq 'seq' || $suffix eq 'clipped'){
	    	my $content = $data_file->slurp;
	    	$content =~ s/^>\S*/> $base/g;

	    	$data_file->spew($content);
	    }
	}
}
