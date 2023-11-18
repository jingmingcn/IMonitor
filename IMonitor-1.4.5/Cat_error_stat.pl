#!/usr/bin/perl -w
use strict;

die "perl $0 <>\n" unless(@ARGV==1);
open I, "$ARGV[0]" or die;
my %all;
while(<I>)
{
	chomp;
	next if(/^$/);
	my @line = split;
	$all{$line[0]} += $line[1];
}
close I;

for("without_V_or_J:","V_and_J_strand_conflict:","CDR3_length_error:","all_raw_seq:","high_quality_seq:","low_quality_correct:","low_quality_filter:","low_quality_unmerge:","new_core_sequence:","PCR_sequence_correct:")
{
	print "$_\t$all{$_}\n";
}
