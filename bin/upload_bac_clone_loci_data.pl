#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use LIMS2::Model;
use Perl6::Slurp;
use Log::Log4perl ':easy';
use Getopt::Long;
use Try::Tiny;
use feature qw( say );

use Smart::Comments;

my $log_level = $WARN;
my $commit;

GetOptions(
    debug   => sub { $log_level = $DEBUG },
    verbose => sub { $log_level = $INFO },
    commit  => \$commit,
) and @ARGV == 1
    or die "Usage: $0 FILE\n";

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );

my $model = LIMS2::Model->new( user => 'lims2' );

my @loci = map { chomp; $_ } slurp( $ARGV[0] );

for my $loci_data ( @loci ) {
    $model->txn_do(
        sub {
            try{
                update_clone_loci( $loci_data );
            }
            catch {
                $model->txn_rollback;
                ERROR( $_ );
            };
            unless ( $commit ) {
                DEBUG( 'Rollback' );
                $model->txn_rollback;
            }
        }
    );

    Log::Log4perl::NDC->remove;
}

sub update_clone_loci {
    my $loci_data = shift;
    my ( $chr, $start, $end, $id ) = split( "\t", $loci_data );
    # bac name ( id ) in data file is not in the same format as we store
    # the bac name in LIMS2
    my $clean_id = clean_id( $id );
    Log::Log4perl::NDC->push( $clean_id );

    my $bac_clone = $model->schema->resultset( 'BacClone' )->find(
        { name => $clean_id, },
        { prefetch => 'loci' }
    );

    unless ( $bac_clone ) {
        WARN( 'Can not find bac_clone, CREATING IT' );
        $bac_clone = $model->schema->resultset( 'BacClone' )->create(
            {
                name           => $clean_id,
                bac_library_id => 'black6',
            }
        );
    }

    INFO( "Working On bac_clone " . $bac_clone->id );

    check_coordinates( $bac_clone, $chr, $start, $end );
    add_new_clone_loci( $bac_clone, $chr, $start, $end );

    return;
}

sub add_new_clone_loci {
    my ( $bac_clone, $chr, $start, $end ) = @_;

    my $new_locus_params = {
        assembly  => 'GRCm38',
        chr_name  => $chr,
        chr_start => $start,
        chr_end   => $end,
    };

    my $validated_locus = $model->check_params( $new_locus_params, $model->pspec_create_bac_clone_locus );

    $bac_clone->create_related(
        loci => {
            assembly_id => $validated_locus->{assembly},
            chr_id      => $model->_chr_id_for( @{$validated_locus}{ 'assembly', 'chr_name' } ),
            chr_start   => $validated_locus->{chr_start},
            chr_end     => $validated_locus->{chr_end}
        }
    );

    INFO( "Added GRMc38 loci for clone" );
}

sub check_coordinates {
    my ( $bac_clone, $chr, $start, $end ) = @_;
    my $current_loci;

    my @loci = $bac_clone->loci;
    my $loci_num = scalar( @loci );
    if ( $loci_num == 1 ) {
        $current_loci = shift @loci;
        LOGDIE( 'Current loci assembly: ' . $current_loci->assembly_id )
            unless $current_loci->assembly_id eq 'NCBIM37'
    }
    elsif ( $loci_num > 1 ) {
        LOGDIE( "We have $loci_num loci for this clone" );
    }
    else {
        INFO( "We have no loci for this clone" );
        return;
    }

    my $orig_start = $current_loci->chr_start;
    my $orig_end   = $current_loci->chr_end;
    my $orig_chr   = $current_loci->chr->name;

    unless ( $orig_chr eq $chr ) {
        WARN( "Original chr: $orig_chr, new chr: $chr" );
        return;
    }

    DEBUG( "Original Start: $orig_start, New Start: $start, diff " . abs( $orig_start - $start ) );
    DEBUG( "Original End: $orig_end, New End: $end, diff " . abs( $orig_end - $end ) );
}

sub clean_id {
    my $id = shift;
    $id =~ s/"//g;
    $id =~ s/^ID\s+//g;

    if ( $id =~ /RPCI/ ) {
        my ( $surp, $num, $ident ) = split( '-', $id );
        return 'RP' . $num . '-' . $ident;
    }
    else {
        return $id;
    }
}

=head1 SYNOPSIS

upload_bac_clone_loci_data.pl[options] bac_clone_loci_data.csv

      --verbose         Print informational messages
      --debug           Print debug messages
      --commit          Commit changes, by default rolls back any changes

=head1 DESCRIPTION

Expects a file, with data in following format:
chr,start,end,id

NOTE: Currently set up only expecting one set of loci info for bac clones, from NCBIM37

=headl TODO

Make it work for conditional designs

=cut
