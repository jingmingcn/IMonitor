#!/usr/bin/perl -w
use strict;

die "perl $0 <>\n" unless(@ARGV==1);
open I, "$ARGV[0]" or die;
while(<I>)
{
	chomp;
	my @line = split /\|/,$_;
	if($line[4] ne "CH1"){
		<I>;
		next;
	}
	print ">$line[1]\n";
	chomp(my $seq = <I>);
	print "$seq\n";
}
close I;
