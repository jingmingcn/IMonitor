#!/usr/bin/perl -w
use strict;

die "perl $0 <single.gz><id.back><out>\n" unless(@ARGV==3);


open O,"|gzip >$ARGV[2]" or die;
open I, "gzip -dc $ARGV[1]|" or die;

my %id_back;
while(<I>)
{
	chomp;
	my @line = split;
	push @{$id_back{$line[3]}} ,$line[1];
}
close I;


my %all;
open I, "gzip -dc $ARGV[0]|" or die;
while(<I>)
{
	chomp;
	my @line = split;
	my $id = (split /:/,$line[30])[0];
	
	for(@{$id_back{$id}}){
		$line[30] = $_;
		my $new_line = join "\t",@line;
		my $id_n = $_;
		$id_n =~s/^cp//;
		$id_n =~s/c:[0-9]*$//;
		$all{$id_n} = $new_line;	
	}
}
close I;

for(sort {$a<=>$b} keys %all)
{
	print O "$all{$_}\n";
}

