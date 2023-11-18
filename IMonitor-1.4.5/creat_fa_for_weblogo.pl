#!/usr/bin/perl -w
use strict;
my $in=shift;
my $div=shift||`less $in|sort -k3,3n|tail -1|cut -f3|awk '{print int(\$1/5000)+1}'`;
chomp $div;
#print "$div\n" and exit;
my %hash;
my %pre_cnt;
open I, "sort -k1,1n $in|" or die "Failed in opening $in !!\n";
while (<I>){
	chomp;
	my @F=split;
	$F[2]=int ($F[2]/$div);
	$pre_cnt{$F[0]}=0 if(!exists $pre_cnt{$F[0]});
	if($F[2]>0){
		for my $line(($pre_cnt{$F[0]}+1)..($pre_cnt{$F[0]}+$F[2])){
			$hash{$line}{$F[0]}=$F[1];
		}
	}
	$pre_cnt{$F[0]}+=$F[2];
}
close I;
my $max=`sort -k1,1n $in|tail -1|cut -f1`;
my $min = `sort -k1,1n $in|head -1|cut -f1`;
chomp $max;
for my $line(sort {$a<=>$b} keys %hash){
	print ">$line\n";
	for my $pos($min..$max){
		${$hash{$line}}{$pos}="-" if(!exists ${$hash{$line}}{$pos});
		print "${$hash{$line}}{$pos}";
	}
	print "\n";
}
