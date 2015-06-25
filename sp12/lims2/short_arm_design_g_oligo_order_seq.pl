#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Getopt::Long;
use Log::Log4perl ':easy';
use LIMS2::Model;
use Pod::Usage;
use Text::CSV;
use feature qw( say );

my $log_level = $WARN;
GetOptions(
    'help'          => sub { pod2usage( -verbose => 1 ) },
    'man'           => sub { pod2usage( -verbose => 2 ) },
) or pod2usage(2);

Log::Log4perl->easy_init( { level => $log_level, layout => '%p %x %m%n' } );

my $model = LIMS2::Model->new( user => 'lims2' );

my $file = $ARGV[0];
LOGDIE( 'Must specify a input file' ) unless $file;

my $input_csv = Text::CSV->new();
open ( my $input_fh, '<', $file ) or die( "Can not open $file " . $! );
$input_csv->column_names( @{ $input_csv->getline( $input_fh ) } );

say 'well,old_design,gene,int_vector,shortened_arm_design_id,G5_oligo,G3_oligo';
while ( my $data = $input_csv->getline_hr( $input_fh ) ) {
    my $design_id = $data->{shortened_arm_design_id};
    my $design = $model->c_retrieve_design( { id => $design_id } );
    my $oligo_order_seqs = $design->oligo_order_seqs;
    #say $data->{well} . ',' . $data->{gene},
          #','
        #. $data->{design_well} . ','
        #. $design_id . ','
        #. $oligo_order_seqs->{G5} . ','
        #. $oligo_order_seqs->{G3};
    say join( ',', ( 
            $data->{well},
            $data->{old_design},
            $data->{gene},
            $data->{int_vector},
            $design_id,
            $oligo_order_seqs->{G5},
            $oligo_order_seqs->{G3},
        )
    ); 
}
close $input_fh;

__END__

=head1 NAME

short_arm_design_g_oligo_order_seq.pl - G5 and G3 oligo order sequences for short arm designs 

=head1 SYNOPSIS

  short_arm_design_g_oligo_order_seq.pl [options] <file>

      --help            Display a brief help message
      --man             Display the manual page

=head1 DESCRIPTION

=cut
