#!/usr/bin/perl -w
use strict;

die "perl $0 <Basic_Statistics_of_Sequencing_Quality.txt> <sample_raw_data_filter_1.txt(RHUMdqyNTANSRAAPEI-87_raw_data_filter_1.txt)>" unless (@ARGV == 2);

open IN, $ARGV[0] or die $!;
open OUT, ">$ARGV[1]" or die $!;

while(<IN>){
	chomp;
	if($_ =~ /Total/ && $_ =~ /reads/){
		my @tmp = split;
		my $perc = $tmp[6]/$tmp[4];
		print OUT "All_read_number(PE=1):\t$tmp[4]\nretained_read_num:\t$tmp[6]\t$perc\n";
	}


}
