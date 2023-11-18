#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
my($ref_V,$ref_J,$output_V,$output_J);
GetOptions(
	"v:s" => \$ref_V,
	"j:s" => \$ref_J,
	"ov:s" => \$output_V,
	"oj:s" => \$output_J
);

open IN_V,"$ref_V" or die "can't open the file $ref_V\n";
open IN_J,"$ref_J" or die "can't open the file $ref_J\n";
open OUT_V,">$output_V" or die "can't open the file $output_V\n";
open OUT_J,">$output_J" or die "can't open the file $output_J\n";
my$count = -1;
while(<IN_J>){
	chomp;
	my@line = split;
	if($line[0] =~ /\>/){
		my$len = length($line[0]);
		my$ref_J_name = (split /\|/,$_)[1];
		print OUT_J "$ref_J_name\t";
	}else{
		my$pre_FR4 = substr($line[0],0,23);
		$pre_FR4 =~ tr/-//ds;
		my$FR4_start = length($pre_FR4) + 1;
		my$FR4 = substr($line[0],23,31);
                $FR4 =~ tr/-//ds;
		$FR4 =~ tr/acgt/ACGT/;	
		my$FR4_end = length($FR4) + length($pre_FR4);
		print OUT_J "$FR4\t$FR4_start\t$FR4_end\n";
	}	
}
close IN_J;
while(<IN_V>){
	chomp;
	my@line = split;
	if($line[0] =~ /\>/){
		my$len = length($line[0]);
		my$ref_name = (split /\|/,$_)[1];
		print OUT_V "$ref_name\t";
	}else{
		my($FR1,$CDR1,$FR2,$CDR2,$FR3);
		$FR1 = substr($line[0],0,78);
		$FR1 =~ tr/.//ds;
		my$FR1_len = length($FR1);
		print OUT_V "$FR1\t1\t$FR1_len\t";

		$CDR1 = substr($line[0],78,36);
		$CDR1 =~ tr/.//ds;
		my$CDR1_len = length($CDR1);
		my$CDR1_start = $FR1_len + 1;
		my$CDR1_end = $FR1_len + $CDR1_len;
		print OUT_V "$CDR1\t$CDR1_start\t$CDR1_end\t";

		$FR2 = substr($line[0],114,51);
		$FR2 =~ tr/.//ds;
		my$FR2_len = length($FR2);
		my$FR2_start = $CDR1_end + 1;
		my$FR2_end = $CDR1_end + $FR2_len;
		print OUT_V "$FR2\t$FR2_start\t$FR2_end\t";

		$CDR2 = substr($line[0],165,30);
		$CDR2 =~ tr/.//ds;
		my$CDR2_len = length($CDR2);
		my$CDR2_start = $FR2_end + 1;
		my$CDR2_end = $FR2_end + $CDR2_len;
		print OUT_V "$CDR2\t$CDR2_start\t$CDR2_end\t";

		$FR3 = substr($line[0],195,117);
		$FR3 =~ tr/.//ds;
		my$FR3_len = length($FR3);
		my$FR3_start = $CDR2_end + 1;
		my$FR3_end = $CDR2_end + $FR3_len;
		print OUT_V "$FR3\t$FR3_start\t$FR3_end\t";
		
		print OUT_V "$FR1$CDR1$FR2$CDR2$FR3\n";
	}
	
}
close IN_V;
