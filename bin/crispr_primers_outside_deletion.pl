#!/usr/bin/env perl

use strict;
use warnings;

use LIMS2::Model;
use Log::Log4perl qw( :easy );
use Data::Dumper;
use feature qw(say);
use Try::Tiny;


=head

Ticket #12371

Some crispr sequencing primers (SF1/SR1) fall outside the genomic region that is deleted by the design.
In these cases the sequencing done for crispr-QC picks up both those with crispr damage and those with
the designed deletion.

Need to identify all assembly wells/experiments where the crispr entity primers fall outside
the design deleted region.

=cut


open (my $fh, ">", "assembly_well_crispr_primer_check.txt") or die $!;
my @result_items = qw(
    primer_name
    primer_id
    outsite_deletion
    primer_start
    primer_end
    target_start
    target_end
    design_id
);

print $fh join "\t", qw(well_name experiment_id), @result_items;
print $fh "\n";

my $model = LIMS2::Model->new({ user => 'lims2' });

my @wells = map { $_->wells }
            $model->schema->resultset('Plate')->search({ type_id => 'ASSEMBLY'})->all;

my %experiment_check_result;
say "Checking crispr primer locations for ".scalar(@wells)." assembly wells";

foreach my $well (@wells){
    foreach my $exp ($well->experiments){
    	my $check_result = $experiment_check_result{ $exp->id };
    	unless ($check_result){
    		# Do the check and add it to hash
    		$check_result = _do_experiment_primer_check($exp);
    		$experiment_check_result{ $exp->id } = $check_result;
    	}

        # Output result
        foreach my $result (@{ $check_result }){
	        my @output_items = ($well->as_string, $exp->id);
	        push @output_items, map { $result->{$_} } @result_items;
	        print $fh join "\t", @output_items;
	        print $fh "\n";
        }
    }
}

sub _do_experiment_primer_check{
	my ($exp) = @_;

    say "checking experiment ".$exp->id;

    my @result;
    my $design = $exp->design;

    # Check is not relevant to crispr only experiments
    return \@result unless $design;

    my @primers = map { $exp->crispr_entity->current_primer($_) } qw(SF1 SR1);

	foreach my $primer (@primers){
		next unless $primer;

		say "checking primer ".$primer->primer_name->primer_name;

        my $start = $primer->start;
        my $end = $primer->end;

		my $primer_result = {
			primer_name => $primer->primer_name->primer_name,
			primer_id   => $primer->id,
			primer_start => $start,
			primer_end   => $end,
			design_id    => $design->id,
			target_start => $design->target_region_start,
			target_end   => $design->target_region_end,
		};

        if($start < $design->target_region_start or $end > $design->target_region_end){
            $primer_result->{outsite_deletion} = 1;
        }
        else{
            $primer_result->{outsite_deletion} = 0;
        }
        push @result, $primer_result;
	}
	return \@result;
}