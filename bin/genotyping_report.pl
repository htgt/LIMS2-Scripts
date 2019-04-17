#!/usr/bin/perl
use strict;
use warnings;
use CAM::PDF;
use Carp;
use Getopt::Long;
use LWP::Simple qw/get getstore/;
use PDF::WebKit;
use Pod::Usage;
use Readonly;
our $VERSION = 0.01;
Readonly::Scalar my $BASE_URL => $ENV{LIMS2_GENOTYPING_BASEURL}
  // 'https://www.sanger.ac.uk/htgt/lims2';
Readonly::Hash my %PAGE_COUNT_HANDLERS => (

    # delete typographic orphan
    3 => sub { my $pdf = shift; $pdf->deletePage( $pdf->numPages ); },

    # 2 is fine
    2 => sub { },
);

if ( exists $ENV{LIMS2_GENOTYPING_WKHTMLTOPDF} ) {
    PDF::WebKit->configure(
        sub {
            $_->wkhtmltopdf( $ENV{LIMS2_GENOTYPING_WKHTMLTOPDF} );
        }
    );
}

sub print_page {
    my ( $plate, $well ) = @_;
    my $page     = "public_reports/well_genotyping_info/$plate/$well";
    my $contents = get("$BASE_URL/$page");

    # replace our version of bootstrap with a more modern one
    # (that doesn't mess up the printing of background colours)
    # this should also be robust with the staging URL rewriting bug
    $contents =~ s{[^"]+bootstrap.min.css}
        {https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css}gxms;

    my $kit = PDF::WebKit->new( \$contents, page_size => 'letter' );
    my $output_file = "${plate}_${well}.pdf";
    $kit->to_file($output_file);
    return $output_file;
}

sub append_pdf {
    my ( $original, $file ) = @_;
    my $pdf       = CAM::PDF->new($file);
    my $num_pages = $pdf->numPages;
    if ( exists $PAGE_COUNT_HANDLERS{$num_pages} ) {
        $PAGE_COUNT_HANDLERS{$num_pages}->($pdf);
    }
    else {
        # 1 page is _probably_ a 404, more is... interesting.
        print {*STDERR} "WARNING: $file has $num_pages page(s)\n";
    }

    # do the append, or return the new file if there's no original
    if ($original) {
        $original->appendPDF($pdf);
        return $original;
    }
    return $pdf;
}

my ( $input_file, $help );
my $output_file = 'genotyping_report.pdf';
GetOptions(
    'input=s'  => \$input_file,
    'output=s' => \$output_file,
    'help|?'   => \$help,
) or croak 'Cannot parse arguments';
pod2usage( -exitval => 0, -verbose => 3 ) if $help;
croak "--input argument is required\n" if not $input_file;
croak "input file $input_file not found\n" if not -e $input_file;

my $pdf;
open my $input, '<', $input_file or croak "Cannot open $input_file: $!\n";
while ( my $name = <$input> ) {
    chomp $name;
    my ( $plate, $well ) = $name =~ m/(.+)_([[:upper:]]\d+)$/xms;
    my $file = print_page( $plate, $well );
    $pdf = append_pdf( $pdf, $file );
}
close $input or croak "Cannot close input file: $!\n";

$pdf->cleanoutput($output_file);

__END__

=head1 NAME

genotyping_report.pl - Prints a list of genotyping reports

=head1 USAGE

genotyping_report -i <input file> [-o <output file>]

=head1 OPTIONS

=over 8

=item B<--input>

Specifies a file containing a list of input wells (in the form <PLATE>_<WELL>) to retrieve genotyping reports for.

=item B<--output>

Specifies an output file. Defaults to F<genotyping_report.pdf>

=item B<--help>

Displays this help text.

=back

=head1 DESCRIPTION

B<This program> reads a list of wells and produces a list of PDFs containing all their genotyping reports (plus a combined summary PDF).

=head1 DEPENDENCIES

=over 8

=item wkhtmltopdf

=back

=head1 ENVIRONMENT

=over 8

=item LIMS2_GENOTYPING_BASEURL

The base URL for the LIMS2 setup. Assumes L<https://www.sanger.ac.uk/htgt/lims2> if not specified.

=item LIMS2_GENOTYPING_WKHTMLTOPDF

This script depends on wkhtmltopdf (L<https://wkhtmltopdf.org/>) to output PDFs. If it's not installed in the PATH give the path here instead.

=back

=cut

