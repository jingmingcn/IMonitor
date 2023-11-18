#!/usr/bin/perl -w 
use strict;

die "perl $0 <in-dir><sample><out><num>\n" unless(@ARGV==4);

my %all;
my $num = $ARGV[3];

for(my $i = 1 ;$i<=$num ; $i++)
{
	open I, "$ARGV[0]/$ARGV[1]_VDJ.alignment.stat_2.$i.txt" or die;
	while(<I>)
	{
		chomp;
		my @line = split;
		$all{$line[0]} += $line[1];
	}
	close I;
}

open O, ">$ARGV[2]" or die;
print O "sequence_num:\t$all{'sequence_num:'}\n";
for ("V_alignment_rate:","D_alignment_rate:","J_alignment_rate:","VJ_alignment_rate:","VDJ_alignment_rate:")
{
	my $r = $all{$_}/$all{'sequence_num:'}*100;
	print O "$_\t$all{$_}\t$r\n";
}
