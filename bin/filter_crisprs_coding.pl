#! /usr/bin/perl

use LIMS2::Model;
use feature qw/ say /;
use strict;
use warnings;
use Try::Tiny;
=head filter_crisprs_coding

This file writes out filtered_crisprs.tsv

NB: Requires an adidtional table in LIMS2 that contains the translated exon coordinates for each exon that has a translated region.

This is used as part fot he first with join. The SQL to create the table is in the sql directory. The populate_gene_coding_meta.pl
script is used to generate a csv file for loading the table in psql, using:

psql> \copy gene_coding_meta from gene_coding_meta.csv with csv header;

=cut

my @exon_array = (qw/
ENSMUSE00000333247
ENSMUSE00000350958
ENSMUSE00000464183
ENSMUSE00000488226
ENSMUSE00000196278
ENSMUSE00000196281
ENSMUSE00000404314
ENSMUSE00000422460
ENSMUSE00000422283
ENSMUSE00000396749
ENSMUSE00000333247
ENSMUSE00000296390
ENSMUSE00000652713
ENSMUSE00000489209
ENSMUSE00000557527
ENSMUSE00000538745
ENSMUSE00000230695
ENSMUSE00000263754 ENSMUSE00000263747
ENSMUSE00000470657
ENSMUSE00000492205
ENSMUSE00000365592
ENSMUSE00000448516
ENSMUSE00000608997
ENSMUSE00000512643
ENSMUSE00001051733
ENSMUSE00000456917 ENSMUSE00000388025 
ENSMUSE00000903869 ENSMUSE00000593974 
ENSMUSE00001052490 ENSMUSE00001057978 
ENSMUSE00000593996 ENSMUSE00000531437 
ENSMUSE00000403934 ENSMUSE00000666578 
ENSMUSE00000555722 ENSMUSE00000703719 
ENSMUSE00001066288 ENSMUSE00001085573 ENSMUSE00001055090 ENSMUSE00000978336
ENSMUSE00001022393 ENSMUSE00000100522 ENSMUSE00000100524 ENSMUSE00000880341
ENSMUSE00000513477 ENSMUSE00000371401 ENSMUSE00000330942 ENSMUSE00000330932
ENSMUSE00001246107 ENSMUSE00001293522 ENSMUSE00000198849 ENSMUSE00001267780
ENSMUSE00000889409 ENSMUSE00000584592 ENSMUSE00000584591 ENSMUSE00000896048
ENSMUSE00000811778 ENSMUSE00001249639 ENSMUSE00000817795 ENSMUSE00000817110
ENSMUSE00000459571 ENSMUSE00000459564 ENSMUSE00000459558 ENSMUSE00000331337
ENSMUSE00000674785 ENSMUSE00000275186 ENSMUSE00000674788 ENSMUSE00000674786
ENSMUSE00001092285 ENSMUSE00001074951 ENSMUSE00000980415 ENSMUSE00000477595
ENSMUSE00000452300 ENSMUSE00000364411 
ENSMUSE00000459879 ENSMUSE00000510070 
ENSMUSE00000900523 ENSMUSE00000486395 
ENSMUSE00001246281 ENSMUSE00001230941 
ENSMUSE00000648567 ENSMUSE00000391241 ENSMUSE00000507342
ENSMUSE00000850412 ENSMUSE00000990178 ENSMUSE00000858679
ENSMUSE00000703726 ENSMUSE00001026465 ENSMUSE00000703724
ENSMUSE00000621818 ENSMUSE00000621817 ENSMUSE00000437793
ENSMUSE00000682630 ENSMUSE00000682629 ENSMUSE00000682627
    ENSMUSE00000213365 ENSMUSE00000350978 ENSMUSE00000213363 ENSMUSE00000357180

    ENSMUSE00000877045

    ENSMUSE00000933442 ENSMUSE00000732266

    ENSMUSE00000493401 ENSMUSE00000640012 ENSMUSE00000506216

    ENSMUSE00000245196 ENSMUSE00000245188 ENSMUSE00000245177 ENSMUSE00000363331

    ENSMUSE00001286347 ENSMUSE00000642334 ENSMUSE00000642333 ENSMUSE00000642332 ENSMUSE00001240662
    ENSMUSE00001216320 ENSMUSE00000642328 ENSMUSE00000642327 ENSMUSE00000642325 ENSMUSE00001253463 ENSMUSE00000746080

    ENSMUSE00000653699 ENSMUSE00000653698 ENSMUSE00000292995 ENSMUSE00000292988 ENSMUSE00000713967 ENSMUSE00000455052
    ENSMUSE00000653700 ENSMUSE00001283266 ENSMUSE00000455041 ENSMUSE00001309548 ENSMUSE00001298290 ENSMUSE00001208340
    ENSMUSE00000721479 ENSMUSE00001233499

    ENSMUSE00000643852 ENSMUSE00000643850 ENSMUSE00000643849 ENSMUSE00000643848 ENSMUSE00000643847 ENSMUSE00000643846
    ENSMUSE00000643845 ENSMUSE00000643844 ENSMUSE00000643842 ENSMUSE00000643841 ENSMUSE00000643839 ENSMUSE00000643838
    ENSMUSE00000643837 ENSMUSE00000643836
/);

my $exon_sql_list = join ',', @exon_array;

my $model = LIMS2::Model->new( { user => 'webapp', audit_user => 'dp10@sanger.ac.uk' } );

my $sql_query = five_query( $exon_sql_list);

my $sql_result_5_prime =  $model->schema->storage->dbh_do(
    sub {
         my ( $storage, $dbh ) = @_;
         my $sth = $dbh->prepare_cached( $sql_query );
         $sth->execute( );
         $sth->fetchall_arrayref({
             'ensembl_exon_id' => 1,
             'marker_symbol' => 1,
             'exon_rank' => 1,
             'crispr_pair_id' => 1,
             'left_diff' => 1,
             'right_diff' => 1,
             'pair_length' => 1,
             'spacer' => 1,
             'left_window' => 1,
             'right_window' => 1,
             'left_crispr_seq' => 1,
             'right_crispr_seq' => 1,
             'exon_start' => 1,
             'exon_length' => 1,
             'L0' => 1,
             'L1' => 1,
             'L2' => 1,
             'R0' => 1,
             'R1' => 1,
             'R2' => 1,
         });
    }
);
say ('SQL query for 5_prime_query brought back ' . @{$sql_result_5_prime} . ' rows.' );

$sql_query = three_query( $exon_sql_list);

my $sql_result_3_prime =  $model->schema->storage->dbh_do(
    sub {
         my ( $storage, $dbh ) = @_;
         my $sth = $dbh->prepare_cached( $sql_query );
         $sth->execute( );
         $sth->fetchall_arrayref({
             'ensembl_exon_id' => 1,
             'marker_symbol' => 1,
             'exon_rank' => 1,
             'crispr_pair_id' => 1,
             'left_diff' => 1,
             'right_diff' => 1,
             'pair_length' => 1,
             'spacer' => 1,
             'left_window' => 1,
             'right_window' => 1,
             'left_crispr_seq' => 1,
             'right_crispr_seq' => 1,
             'exon_end' => 1,
             'exon_length' => 1,
             'L0' => 1,
             'L1' => 1,
             'L2' => 1,
             'R0' => 1,
             'R1' => 1,
             'R2' => 1,
         });
    }
);
say ('SQL query for 3_prime_query brought back ' . @{$sql_result_3_prime} . ' rows.' );

# Holding area
my $clip;
# Process 5 prime filtered crisprs into clip
my @headers_5_prime =  (qw/
            ensembl_exon_id
            marker_symbol
            exon_rank
            crispr_pair_id
            left_diff
            right_diff
            pair_length
            left_window
            right_window
            left_crispr_seq
            right_crispr_seq
            exon_start
            exon_length
            spacer
            L0
            L1
            L2
            R0
            R1
            R2
            /);
my @headers_3_prime = (qw/
            ensembl_exon_id
            marker_symbol
            exon_rank
            crispr_pair_id
            left_diff
            right_diff
            pair_length
            left_window
            right_window
            left_crispr_seq
            right_crispr_seq
            exon_end
            exon_length
            spacer
            L0
            L1
            L2
            R0
            R1
            R2
            /);
my $headers_5_tsv = join "\t", @headers_5_prime;
my $headers_3_tsv = join "\t", @headers_3_prime;

foreach my $row ( @{$sql_result_5_prime} ) {
    if ( ! defined $row->{'exon_rank'} ) {
        say $row->{'marker_symbol'} . ':' . $row->{'ensembl_exon_id'} . ' has no rank information - assigned to rank none';
        $row->{'exon_rank'} = 'none';
    }
    my $out_row = join "\t", @{$row}{@headers_5_prime};
    push @{$clip->{$row->{'marker_symbol'}}->{$row->{'exon_rank'}}->{'5_prime'}}, $out_row; 
}
foreach my $row ( @{$sql_result_3_prime} ) {
    if ( ! defined $row->{'exon_rank'} ) {
        say $row->{'marker_symbol'} . ':' . $row->{'ensembl_exon_id'} . ' has no rank information - assigned to rank none';
        $row->{'exon_rank'} = 'none';
    }
    my $out_row = join "\t", @{$row}{@headers_3_prime};
    push @{$clip->{$row->{'marker_symbol'}}->{$row->{'exon_rank'}}->{'3_prime'}}, $out_row; 
}
my $tab_filename = 'filtered_crisprs.tsv';
open( my $tab_fh, '>', $tab_filename )
    or die "Can't open file $tab_filename: $! \n";

foreach my $gene_symbol ( keys %{$clip} ) {
    print $tab_fh "$gene_symbol: 5_prime\n"; 
    print $tab_fh $headers_5_tsv . "\n";
    foreach my $exon_rank ( sort keys %{$clip->{$gene_symbol}} ) {
        if ( $clip->{$gene_symbol}->{$exon_rank}->{'5_prime'} ) {
            print $tab_fh "\n";
            foreach my $line ( @{$clip->{$gene_symbol}->{$exon_rank}->{'5_prime'}} ) {
                if ( defined($line) and length($line) ) {
                    print $tab_fh $line . "\n";
                }
                else
                {
                    print $tab_fh 'No results with selected parameters' . "\n";
                }
            }
        }
        else
        {
            print $tab_fh 'No results with selected parameters for this exon' . "\n";
        }
    }
    print $tab_fh "\n";
    print $tab_fh "$gene_symbol: 3_prime\n"; 
    print $tab_fh $headers_3_tsv . "\n";
    
    foreach my $exon_rank ( sort keys %{$clip->{$gene_symbol}} ) {
        if ( $clip->{$gene_symbol}->{$exon_rank}->{'3_prime'} ) {
            print $tab_fh "\n";
            foreach my $line ( @{$clip->{$gene_symbol}->{$exon_rank}->{'3_prime'}} ) {
                if ( defined($line) and length($line) ) {
                print $tab_fh $line . "\n"
                }
                else
                {
                    print $tab_fh 'No results with selected parameters' . "\n";
                }
            }
        }
        else
        {
            print $tab_fh 'No results with selected parameters for this exon' . "\n";
        }
}
    print $tab_fh "\n\n";
}
print $tab_fh "End of File\n";
close $tab_fh;

say 'Done' . "\n";
##
sub five_query {
    my  $exon_sql_list = shift;

return <<"EOQ";

WITH dt as (
    SELECT
         des.ensembl_exon_id
        ,des.marker_symbol
        ,exon_rank
        ,chr_id
        ,(gc.cd_start-200) as chr_start
        ,gc.cd_start as exon_start
        ,(gc.cd_end+200) as chr_end
        ,cd_end as exon_end
        ,gc.cd_end-gc.cd_start as exon_length
        ,100 as left_window
        ,100 as right_window
    FROM (SELECT unnest('
    {$exon_sql_list}
    '::text[]) AS id) x
    JOIN design_targets des ON des.ensembl_exon_id=x.id
    JOIN gene_coding_meta gc ON des.ensembl_exon_id=gc.ensembl_exon_id
),
pairs as (
select b.id "pair_id"
	, a.crispr_id "left_crispr_id"
	, d.seq "left_crispr_seq"
	, a.chr_start "left_crispr_start"
	, a.chr_end "left_crispr_end"
	, c.crispr_id "right_crispr_id"
	, e.seq "right_crispr_seq"
	, c.chr_start "right_crispr_start"
	, c.chr_end "right_crispr_end"
	, b.spacer
	
	from dt
	join crispr_loci a on a.chr_start >= dt.chr_start and a.chr_end <= dt.chr_end
        and a.chr_id = dt.chr_id
	join crispr_pairs b on a.crispr_id = b.left_crispr_id
	join crispr_loci c on b.right_crispr_id = c.crispr_id
	join crisprs d on b.left_crispr_id = d.id
	join crisprs e on b.right_crispr_id = e.id
),
exon_left as (
SELECT
     dt.ensembl_exon_id
    ,dt.marker_symbol
    ,dt.exon_rank
    ,dt.left_window
    ,dt.right_window
    ,pairs.pair_id crispr_pair_id
    ,pairs.left_crispr_id
    ,pairs.right_crispr_id
    ,dt.exon_start exon_start
    ,dt.exon_length
    ,pairs.left_crispr_start-dt.exon_start as lleft_diff
    ,pairs.right_crispr_end-dt.exon_start as lright_diff
    ,pairs.left_crispr_start
    ,pairs.left_crispr_end
    ,pairs.right_crispr_start
    ,pairs.right_crispr_end
    ,pairs.left_crispr_seq
    ,pairs.right_crispr_seq
    ,pairs.spacer
    FROM dt
    JOIN pairs
    ON left_crispr_start > dt.exon_start - dt.left_window
        AND right_crispr_end < dt.exon_start + dt.right_window
),
results as (
SELECT
     ensembl_exon_id
    ,marker_symbol
    ,exon_rank
    ,crispr_pair_id
    ,lleft_diff as "left_diff"
    ,lright_diff as "right_diff"
    ,abs(lleft_diff) + abs(lright_diff) as pair_length
    ,spacer
    ,left_window
    ,right_window
    ,left_crispr_seq
    ,right_crispr_seq
    ,exon_start
    ,exon_length
    ,split_part(split_part(coff_left.summary,',',1), ':', 2)::Integer as L0
    ,split_part(split_part(coff_left.summary,',',2), ':', 2)::Integer as L1
    ,split_part(split_part(coff_left.summary,',',3), ':', 2)::Integer as L2
    ,split_part(split_part(coff_right.summary,',',1), ':', 2)::Integer as R0
    ,split_part(split_part(coff_right.summary,',',2), ':', 2)::Integer as R1
    ,split_part(split_part(coff_right.summary,',',3), ':', 2)::Integer as R2
FROM exon_left

JOIN crispr_off_target_summaries coff_left
    ON coff_left.crispr_id=exon_left.left_crispr_id

JOIN crispr_off_target_summaries coff_right
    ON coff_right.crispr_id=exon_left.right_crispr_id
)

SELECT DISTINCT * FROM results
WHERE l0 = 1 AND L1 = 0 and L2 <= 1
    AND R0 = 1 AND R1 = 0 and R2 <= 1
    AND spacer >= 10 and spacer <= 20
    ORDER BY results.marker_symbol, results.exon_rank, spacer
EOQ
}

sub three_query {
    my  $exon_sql_list = shift;

return <<"EOQ";
WITH dt as (
    SELECT
         des.ensembl_exon_id
        ,des.marker_symbol
        ,exon_rank
        ,chr_id
        ,(gc.cd_start-200) as chr_start
        ,gc.cd_start as exon_start
        ,(gc.cd_end+200) as chr_end
        ,cd_end as exon_end
        ,gc.cd_end-gc.cd_start as exon_length
        ,100 as left_window
        ,100 as right_window
    FROM (SELECT unnest('
    {$exon_sql_list}
    '::text[]) AS id) x
    JOIN design_targets des ON des.ensembl_exon_id=x.id
    JOIN gene_coding_meta gc ON des.ensembl_exon_id=gc.ensembl_exon_id
),
pairs as (
select b.id "pair_id"
	, a.crispr_id "left_crispr_id"
	, d.seq "left_crispr_seq"
	, a.chr_start "left_crispr_start"
	, a.chr_end "left_crispr_end"
	, c.crispr_id "right_crispr_id"
	, e.seq "right_crispr_seq"
	, c.chr_start "right_crispr_start"
	, c.chr_end "right_crispr_end"
	, b.spacer
	
	from dt
	join crispr_loci a on a.chr_start >= dt.chr_start and a.chr_end <= dt.chr_end
        and a.chr_id = dt.chr_id
	join crispr_pairs b on a.crispr_id = b.left_crispr_id
	join crispr_loci c on b.right_crispr_id = c.crispr_id
	join crisprs d on b.left_crispr_id = d.id
	join crisprs e on b.right_crispr_id = e.id
),
exon_right as (
SELECT
     dt.ensembl_exon_id
    ,dt.marker_symbol
    ,dt.exon_rank
    ,dt.left_window
    ,dt.right_window
    ,pairs.pair_id crispr_pair_id
    ,pairs.left_crispr_id
    ,pairs.right_crispr_id
    ,dt.exon_end exon_end
    ,dt.exon_length
    ,pairs.left_crispr_start-dt.exon_end as rleft_diff
    ,pairs.right_crispr_end-dt.exon_end as rright_diff
    ,dt.exon_start
    ,pairs.left_crispr_start
    ,pairs.left_crispr_end
    ,pairs.right_crispr_start
    ,pairs.right_crispr_end
    ,pairs.left_crispr_seq
    ,pairs.right_crispr_seq
    ,pairs.spacer
    FROM dt
    JOIN pairs
    ON left_crispr_start > dt.exon_end - dt.left_window
        AND right_crispr_end < dt.exon_end + dt.right_window
),
results as (
SELECT
     ensembl_exon_id
    ,marker_symbol
    ,exon_rank
    ,crispr_pair_id
    ,rleft_diff as "left_diff"
    ,rright_diff as "right_diff"
    ,abs(rleft_diff) + abs(rright_diff) as pair_length
    ,spacer
    ,left_window
    ,right_window
    ,left_crispr_seq
    ,right_crispr_seq
    ,exon_end
    ,exon_length
    ,split_part(split_part(coff_left.summary,',',1), ':', 2)::Integer as L0
    ,split_part(split_part(coff_left.summary,',',2), ':', 2)::Integer as L1
    ,split_part(split_part(coff_left.summary,',',3), ':', 2)::Integer as L2
    ,split_part(split_part(coff_right.summary,',',1), ':', 2)::Integer as R0
    ,split_part(split_part(coff_right.summary,',',2), ':', 2)::Integer as R1
    ,split_part(split_part(coff_right.summary,',',3), ':', 2)::Integer as R2
FROM exon_right

JOIN crispr_off_target_summaries coff_left
    ON coff_left.crispr_id=exon_right.left_crispr_id

JOIN crispr_off_target_summaries coff_right
    ON coff_right.crispr_id=exon_right.right_crispr_id
)
SELECT DISTINCT * from results
WHERE l0 = 1 AND L1 = 0 and L2 <= 1
    AND R0 = 1 AND R1 = 0 and R2 <= 1
    AND spacer >= 10 and spacer <= 20
    ORDER BY results.marker_symbol, results.exon_rank, spacer
EOQ
}


=head old runs
my @exon_array = (qw/
ENSMUSE00000333247
ENSMUSE00000350958
ENSMUSE00000464183
ENSMUSE00000488226
ENSMUSE00000196278
ENSMUSE00000196281
ENSMUSE00000404314
ENSMUSE00000422460
ENSMUSE00000422283
ENSMUSE00000396749
ENSMUSE00000333247
ENSMUSE00000296390
ENSMUSE00000652713
ENSMUSE00000489209
ENSMUSE00000557527
ENSMUSE00000538745
ENSMUSE00000230695
ENSMUSE00000263754 ENSMUSE00000263747
ENSMUSE00000470657
ENSMUSE00000492205
ENSMUSE00000365592
ENSMUSE00000448516
ENSMUSE00000608997
ENSMUSE00000512643
ENSMUSE00001051733
ENSMUSE00000456917 ENSMUSE00000388025 
ENSMUSE00000903869 ENSMUSE00000593974 
ENSMUSE00001052490 ENSMUSE00001057978 
ENSMUSE00000593996 ENSMUSE00000531437 
ENSMUSE00000403934 ENSMUSE00000666578 
ENSMUSE00000555722 ENSMUSE00000703719 
ENSMUSE00001066288 ENSMUSE00001085573 ENSMUSE00001055090 ENSMUSE00000978336
ENSMUSE00001022393 ENSMUSE00000100522 ENSMUSE00000100524 ENSMUSE00000880341
ENSMUSE00000513477 ENSMUSE00000371401 ENSMUSE00000330942 ENSMUSE00000330932
ENSMUSE00001246107 ENSMUSE00001293522 ENSMUSE00000198849 ENSMUSE00001267780
ENSMUSE00000889409 ENSMUSE00000584592 ENSMUSE00000584591 ENSMUSE00000896048
ENSMUSE00000811778 ENSMUSE00001249639 ENSMUSE00000817795 ENSMUSE00000817110
ENSMUSE00000459571 ENSMUSE00000459564 ENSMUSE00000459558 ENSMUSE00000331337
ENSMUSE00000674785 ENSMUSE00000275186 ENSMUSE00000674788 ENSMUSE00000674786
ENSMUSE00001092285 ENSMUSE00001074951 ENSMUSE00000980415 ENSMUSE00000477595
ENSMUSE00000452300 ENSMUSE00000364411 
ENSMUSE00000459879 ENSMUSE00000510070 
ENSMUSE00000900523 ENSMUSE00000486395 
ENSMUSE00001246281 ENSMUSE00001230941 
ENSMUSE00000648567 ENSMUSE00000391241 ENSMUSE00000507342
ENSMUSE00000850412 ENSMUSE00000990178 ENSMUSE00000858679
ENSMUSE00000703726 ENSMUSE00001026465 ENSMUSE00000703724
ENSMUSE00000621818 ENSMUSE00000621817 ENSMUSE00000437793
ENSMUSE00000682630 ENSMUSE00000682629 ENSMUSE00000682627
/);
=cut
