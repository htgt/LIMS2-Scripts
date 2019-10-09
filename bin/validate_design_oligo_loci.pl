#!/usr/bin/perl
use strict;
use warnings;

use Data::Dumper;
use LIMS2::Model;
use LIMS2::Util::EnsEMBL;

sub rename_loci_attr {
    my ($loci, $rename, $original) = @_;

    $loci->{ $rename->{$original} } = $loci->{$original};
    delete $loci->{$original};

    return $loci;
}

my $model  = LIMS2::Model->new( { user => 'lims2' } );

my $ensembl_util = LIMS2::Util::EnsEMBL->new( { 'species' => 'Human' } );
my $slice_adaptor = $ensembl_util->slice_adaptor( 'Human' );

my $design_rs = $model->schema->resultset('Design')->search({
    'me.design_type_id' => 'miseq-nhej'
}, {
    prefetch => {
        'oligos' => 'loci',
    }
});

my $locus;
my @abnormal_oligo_ids;
my @correct;

while (my $design_focus = $design_rs->next) {
    my $oligos = $design_focus->oligos;

    while (my $oligo = $oligos->next) {
        my $loci = $oligo->loci->first;
        $locus->{$design_focus->id}->{$oligo->design_oligo_type_id} = {
            loci => $loci->as_hash,
            seq => {
                db => $oligo->seq,
            },
        };
        $locus->{$design_focus->id}->{$oligo->design_oligo_type_id}->{loci}->{design_oligo_id} = $loci->design_oligo_id;
        $locus->{$design_focus->id}->{$oligo->design_oligo_type_id}->{loci}->{chr_id} = $loci->chr_id;
    }
}

my @design_ids = keys %{ $locus };
foreach my $id (@design_ids) {
    foreach my $name (keys %{ $locus->{$id} }) {
        my $loci = $locus->{$id}->{$name}->{loci};
        print Dumper "$id:$name";
        if ($id == 1016386 && $name eq 'INR'){ 
            $DB::single=1;
        }
        my $slice = $slice_adaptor->fetch_by_region(
            'chromosome',
            $loci->{chr_name},
            $loci->{chr_start},
            $loci->{chr_end},
            $loci->{chr_strand},
        );

        my $ensembl_lookup = $slice->seq;
        my $db_seq = $locus->{$id}->{$name}->{seq}->{db};

        my $oligo_id = "$id:$name";
        if ($ensembl_lookup ne $db_seq) {
            if (length($db_seq) >= length($ensembl_lookup)) {
                push (@abnormal_oligo_ids, {
                    id => $oligo_id,
                    reason => "DB Seq gt or eq. Match not found",
                });
            }
            my $right_trim = $ensembl_lookup;
            chop $right_trim;
            if ($db_seq eq $right_trim) {
                my $amended_end = $loci->{chr_end} - 1;
                my $update = $loci;
                $update->{chr_end} = $amended_end;
                my %rename = (
                    assembly    => 'assembly_id',
                );
                foreach my $key (keys %rename) {
                    $update = rename_loci_attr($update, \%rename, 'assembly');
                }
                delete $update->{species};
                delete $update->{chr_name};
                $model->update_design_oligo_loci($update);
            }
        }
        else {
            push @correct, $oligo_id;
        }
    }
}


print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Abnormal oligos\n";
print Dumper @abnormal_oligo_ids;
print Dumper "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Correct\n";
print Dumper @correct;

1;