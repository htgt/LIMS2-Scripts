#! /usr/bin/perl

use LIMS2::Model;
use feature qw/ say /;
use strict;
use warnings;
use Try::Tiny;

use Bio::EnsEMBL::Registry;

=head Populate_gene_coding_meta

For a given gene with ensembl stable id, return the coding regions from ensembl and store them
in the gene_coding_meta table in LIMS2

=cut

my @gene_array = (qw/
    Gpr27
Rxfp3
Gpr15
Gpr25
Npbwr1
Gpr12
Tas2r131
Tas2r120
Tas2r126
Tas2r135
Tas2r121
Tas2r102
Taar8a
Mrgprb2
Mrgprb1
Mrgpra1
Mrgprx1
Gm4922
Gm7168
Nim1
Npffr1
Npffr2
Sbk2
Scn2b
Gm1078
Gpr82
Kcnc3
Kcng2
Gpr39
Gpr61
Kcns2
Xcr1
Gjc3
Smok2b
Smok4a
Gpr173
Gpr68
Nr2f6
Gpr62
Pak6
Cdc42bpg
Dyrk2
A530099J19Rik
Kcns1
Cdkl3
/);

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
        -host => $ENV{LIMS2_ENSEMBL_HOST} || 'ensembldb.internal.sanger.ac.uk',
        -user => $ENV{LIMS2_ENSEMBL_USER} || 'anonymous'
    );

my $slice_adaptor = $registry->get_adaptor('Mouse', 'Core', 'Slice');

my $model = LIMS2::Model->new( { user => 'webapp', audit_user => 'dp10@sanger.ac.uk' } );

my $gene_clip;

foreach my $gene_symbol ( @gene_array ) {
    $gene_clip->{$gene_symbol} = get_coding_data_canonical( $model, $gene_symbol )    
}
my @headers_csv = (qw/
            marker_symbol
            ensembl_exon_id
            cd_start
            cd_end
            exon_start
            exon_end
            strand
            /);
my $header_string = join ',' , @headers_csv;
my $csv_filename = 'gene_coding_meta_loader.csv';
open( my $csv_fh, '>', $csv_filename )
    or die "Can't open file $csv_filename: $! \n";

print $csv_fh $header_string . "\n";

foreach my $gene_symbol ( keys %{$gene_clip} ) {
    foreach my $exon_id (keys %{$gene_clip->{$gene_symbol}}) {
        my @t_out;
        push @t_out, $gene_symbol;
        push @t_out, $exon_id;
        push @t_out, $gene_clip->{$gene_symbol}->{$exon_id}->{'cd_start'};
        push @t_out, $gene_clip->{$gene_symbol}->{$exon_id}->{'cd_end'};
        push @t_out, $gene_clip->{$gene_symbol}->{$exon_id}->{'exon_start'};
        push @t_out, $gene_clip->{$gene_symbol}->{$exon_id}->{'exon_end'};
        push @t_out, $gene_clip->{$gene_symbol}->{$exon_id}->{'strand'};
        print $csv_fh (join ',', @t_out) . "\n"; 
    }
}


say 'Done..';
exit();

sub get_coding_data_canonical {
    my $model = shift;
    my $gene_symbol = shift; 

    my $rs = $model->schema->resultset('DesignTarget')->search ( { 'marker_symbol' => $gene_symbol } );
    my $ensembl_gene_id = $rs->first->ensembl_gene_id;

    my $ret_val;
    my $slice = $slice_adaptor->fetch_by_gene_stable_id( $ensembl_gene_id, 0 );
    my $genes = $slice->get_all_Genes;
    for my $gene ( @{$genes} ) {
        next if $gene->display_id ne $ensembl_gene_id;
        say 'Gene : ' . $gene->display_id;
        my $transcript = $gene->canonical_transcript;
        say 'Canonical Transcript: ' . $transcript->display_id;
        my $transformed = $transcript->transform('chromosome');
        die if ( ! defined $transformed );
        foreach my $exon ( @{$transcript->get_all_Exons} ) {

            my $cd_start = $exon->coding_region_start( $transcript ); 
            my $cd_end = $exon->coding_region_end( $transcript );
            my $exon_rs = $model->schema->resultset('DesignTarget')->search( { 'ensembl_exon_id' => $exon->display_id } );
            my $lims2_exon_start = $exon_rs->first->chr_start;
            my $lims2_exon_end  = $exon_rs->first->chr_end;
            my $lims2_exon_strand = $exon_rs->first->chr_strand;
            if (defined $cd_start ) {
                say 'Found coding exon: ' . $exon->display_id .  ' starting at: ' . $cd_start . ' ending at: ' . $cd_end;
                $ret_val->{$exon->display_id}->{cd_start} = $cd_start;
                $ret_val->{$exon->display_id}->{cd_end} = $cd_end;
                $ret_val->{$exon->display_id}->{exon_start} = $lims2_exon_start;
                $ret_val->{$exon->display_id}->{exon_end} = $lims2_exon_end;
                $ret_val->{$exon->display_id}->{strand} = $lims2_exon_strand;
                
            }
            else {
                say 'Non-coding exon found';
            }
        }
    }
    return $ret_val;
}

=head
sub get_coding_data_all {
my $slice = $slice_adaptor->fetch_by_gene_stable_id( 'ENSMUSG00000059852', 0 );

my $genes = $slice->get_all_Genes;

for my $gene ( @{$genes} ) {
    say 'Gene : ' . $gene->display_id;
    my $transcripts = $gene->get_all_Transcripts;
    foreach my $transcript ( @{$transcripts} ) {
        say 'Transcript: ' . $transcript->display_id;
        foreach my $exon ( @{$transcript->get_all_Exons} ) {

            my $cd_start = $exon->coding_region_start( $transcript ); 
            if (defined $cd_start ) {
                say 'Found coding exon: ' . $exon->display_id .  ' starting at: ' . $cd_start;
            }
            else {
                say 'Non-coding exon found';
            }
        }
    }
}
return;
}
=cut 


