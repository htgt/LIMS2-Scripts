#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Getopt::Long;
use Pod::Usage;
use Try::Tiny;
use Path::Class;
use feature qw( say );

die "Usage: $0 <species> <plate_id>" unless @ARGV == 2;

my $model = LIMS2::Model->new( user => 'webapp', audit_user => $ENV{USER}.'@sanger.ac.uk' );

{
    my $species = lc($ARGV[0]);
    my $plate_stem;
    if ( $species eq 'mouse' ) {
        $plate_stem = 'M';
    }
    elsif ( $species eq 'human') {
        $plate_stem = 'H';
    }
    else {
        die "Invalid species: " . $ARGV[0];
    }
    my $plate_id = $ARGV[1];
    #cache filehandles. we will need lr OR crisprs, so this way we only make one when needed
    my %filehandles;
    my %filenames = (
        lcrispr  => $plate_stem . "CL" . $plate_id . ".csv",
        rcrispr  => $plate_stem . "CR" . $plate_id . ".csv",
        design   => $plate_stem . "G" . $plate_id . ".csv",
        crispr   => $plate_stem . "C" . $plate_id . ".csv",
        #lcrispr => 'left_crisprs.csv',
        #rcrispr => 'right_crisprs.csv',
        #crispr  => 'crisprs.csv',
        #design  => 'designs.csv',
    );

    sub get_fh {
        my $file = shift;

        die "$file doesnt exist in filenames hash"
            unless exists $filenames{ $file };

        my $filename = $filenames{ $file };

        unless ( defined $filehandles{ $filename } ) {
            say STDERR "Opening file " . $filename;
            $filehandles{ $filename } = file( $filename )->openw;
            print { $filehandles{ $filename } } "well_name,";

            if ( $file =~ /crispr/ ) {
                print { $filehandles{ $filename } } "crispr_id\n";
            }
            elsif ( $file =~ /design/ ) {
                print { $filehandles{ $filename } } "design_id\n";
            }
            else {
                die "Can't identify header name for $filename";
            }
        }

        return $filehandles{ $filename };
    }
}

sub create_wells {
    my @wells; 
    for my $letter ( 'A' .. 'H' ) {
        #use sprintf to pad numbers less than 10 with a 0 at the front
        push @wells, map { $letter . sprintf( "%02d", $_ ) } ( 1 .. 12 ); 
    }

    return @wells;
} 

my @wells = create_wells();
my @crispr_designs = $model->schema->resultset("CrisprDesign")->search(
    { plated => 0 },
    { prefetch => "crispr_pair" }
);

say "Found " . scalar( @crispr_designs ) . " crispr designs";

die "More CrisprDesigns than wells!" if @wells < @crispr_designs;

my $count = 0;
for my $well ( @wells ) {
    if ( ++$count > 48 ) {
        say STDERR "Warning: more than 48 entries in crispr_design table, skipping";
        last;
    }
    my $cd = shift @crispr_designs;

    say { get_fh( 'design' ) } join ",", $well, $cd->design_id;

    if ( defined $cd->crispr_pair_id ) {
        say { get_fh( 'lcrispr' ) } join ",", $well, $cd->crispr_pair->left_crispr_id; 
        say { get_fh( 'rcrispr' ) } join ",", $well, $cd->crispr_pair->right_crispr_id;
    }
    elsif ( defined $cd->crispr_id ) {
        say { get_fh( 'crispr' ) } join ",", $well, $cd->crispr_id;
    }
    elsif ( defined $cd->crispr_group_id ) {
        my @left  = @{ $cd->crispr_group->left_ranked_crisprs };
        my @right = @{ $cd->crispr_group->right_ranked_crisprs };

        die "More than one left crispr in group, code needs changing" unless @left == 1;
        die "More than one right crispr in group, code needs changing" unless @right == 1;

        say { get_fh( 'lcrispr' ) } join ",", $well, $left[0]->id;
        say { get_fh( 'rcrispr' ) } join ",", $well, $right[0]->id;
    }
    else {
        die "Crispr design " . $cd->id . " has no crispr, pair or group linked to it.";
    }

    $cd->update( { plated => 1 } );

    #if we ran out of designs then stop
    last unless @crispr_designs;
}
