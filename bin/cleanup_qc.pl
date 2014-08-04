#!/usr/bin/env perl

use strict;
use warnings;

use Path::Class;
use Getopt::Long;
use Pod::Usage;
use IPC::System::Simple qw( system );

#as team87 run:
#perl /nfs/users/nfs_a/ah19/cleanup_qc.pl --dir qc --do-it
#perl /nfs/users/nfs_a/ah19/cleanup_qc.pl --dir lims2_qc --do-it

#make sure they specified SOMETHING
pod2usage(2) unless @ARGV;

my $user_dir;
my $help = 0;
my $delete = 0; #dry run by default
my $max_per_folder = 50;

GetOptions(
    "help|?"           => \$help,
    "dir=s"            => \$user_dir,
    "do-it"            => \$delete,
    "max-per-folder=i" => \$max_per_folder,
) || pod2usage(2);
pod2usage(1) if $help;

#only allow deletion from directories specified in this hash
my %valid_dirs = (
    qc            => 'qc',
    prescreen_qc  => 'prescreen_qc',
    lims2_qc      => 'lims2_qc',
    lims2_designs => 'lims2_designs'
);

#make sure we only get a valid folder, wouldn't want someone accidentally deleting loads of stuff
die "The directory to clean must be one of: " . join( ", ", keys %valid_dirs ) . "\n"
    unless exists $valid_dirs{ $user_dir };

my $dir = dir("/", "lustre", "scratch110", "sanger", "team87", $valid_dirs{ $user_dir } );
my @runs = get_all_runs( $dir );

#make sure there are enough runs in the folder to delete anything
die "Less than 50 runs present: " . scalar @runs . ". Nothing will be deleted.\n"
    unless @runs > $max_per_folder;

#remove all runs we want to keep from the list (the first runs are the oldest) 
my @to_delete = @runs[ 0 .. (@runs - $max_per_folder) - 1 ];

print scalar @to_delete . " to be deleted from $dir\n";

#loop through the runs, and delete them if --do-it is set.
for my $run ( @to_delete ) {
    my $cmd = sprintf( "rm -rf '%s'", $dir->subdir( $run->{ run_id } ) );
    
    if ( $delete ) {
        print "Deleting " . $run->{ run_id } . "\n";
        system( $cmd );
    }
    else {
        print "Dry run: $cmd\n"
    }
}

#return all sub directories that match our criteria sorted by creation time
sub get_all_runs {
    my ( $dir ) = @_;

    #only add directories that are (hopefully) uuids with a params.yaml file inside.
    #(uuids are always 36 chars long)
    my @child_dirs = grep { $_->is_dir 
                        and -e $_->file( 'params.yaml' )
                        and length( $_->basename ) == 36 } $dir->children;

    #get an array of hashrefs sorted by created time with just the information we want.
    my @runs = sort { $a->{ ctime } <=> $b->{ ctime } }
                    map { 
                            { 
                                run_id => $_->basename, 
                                ctime  => $_->file("params.yaml")->stat->ctime 
                            } 
                    } @child_dirs;

    return @runs;
}

__END__

=head1 NAME

cleanup_qc - delete excess runs from a qc folder

=head1 SYNOPSIS

cleanup_qc.pl --dir <qc|lims2_qc|lims2_designs> [--do-it] [--max-per-folder 50]

Options:

    -dir            The name of the folder you want to clean (must be one of qc|lims2_qc|lims2_designs)
    -do-it          Actually delete the folders, without this it will be a dry run
    -max-per-folder Specify the number of runs to keep within the folder
    -help           Display this message

Example:
    cleanup_qc.pl --dir qc --do-it
    cleanup_qc.pl --dir lims2_qc --max-per-folder 30 --do-it

=head1 DESCRIPTION

B<This program> will loop through all runs within the specified qc folder and delete
any excess folders (determined by --max-per-folder which defaults to 50.)

=cut
