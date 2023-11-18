#!/usr/bin/perl -w
use strict;

=head1 Usage
	perl Merge_filter_fqtofa.pl  <in.fa.gz><out.fa.gz> <change_id.back.gz> <insert-size.distribution>
	<out.fa.gz>/<change_id.back.gz>: must be the suffix ".gz"

=head1 Description
	1. combined the same sequences
	2. insert size length distribution statistic

=cut

die `pod2text $0` unless(@ARGV==4);


my ($inf, $out , $out_id , $out_insert) = @ARGV;

if($inf=~/\.gz$/){
	open I, "gzip -dc $inf|" or die;
}else{
	open I, "$inf" or die;
}

open OUT, "|gzip >$out" or die;
open OD, "|gzip >$out_id" or die;
#open OF, "|gzip >$out_filter" or die;

my %insertsize_len;
#my $java_meger_n = 0;
#my $filter_n = 0;
my $flag = 0;
my $sum = 0;
my $flag_c = 0;
my %combine_seq;
my $total = 0;

# 1. read the cope-mereged file
while(<I>)
{
	chomp;
	my $id = $_;
	chomp(my $seq = <I>);
	$total++;
	# 1.1 insert size length 
	my $len = length($seq);
	$insertsize_len{$len}++;
	$sum++;


	# 1.3 change fq to fa
	$flag++;
	my $new_id = "cp$flag";
	$id = (split)[0];
	my $com_id;
	if(exists $combine_seq{$seq}){
		my ($d , $n) = split /:/,$combine_seq{$seq};
		$n++;
		$combine_seq{$seq} = "$d:$n";
		$com_id = $d;
	}else{
		$flag_c++;
		$com_id = "cp${flag_c}c";
		$combine_seq{$seq} = "$com_id:1";
	}
	print OD "$id\t$new_id\t$len\t$com_id\n";
#	print OUT ">$new_id\n$seq\n";
}
close I;



# output the combined sequence
for(keys %combine_seq)
{
	print OUT ">$combine_seq{$_}\n$_\n";
}

# output insert size length distribution
open O, ">$out_insert" or die;
for(sort {$a<=>$b} keys %insertsize_len)
{
	my $rate = $insertsize_len{$_}/$sum*100;
	print O "$_\t$insertsize_len{$_}\t$rate\n";
}
close O;

#my $tmp = `awk 'NR==3{print \$1"\t"\$2"\t"\$3}' $log`;
#chomp($tmp);
#my ($total, $cope_n , $cope_r) = split /\t/,$tmp;
#my $unmerge_n = $total-$cope_n-$java_meger_n;

print "Total_seq: $total\n";
#print "Merged_seq_raw: ",$cope_n+$java_meger_n,"\t",($cope_n+$java_meger_n)/$total*100,"\n";
#print "\tCOPE_merged_seq_raw: $cope_n\t",$cope_n/$total*100,"\n";
#print "\tJava_merged_seq_raw: $java_meger_n\t",$java_meger_n/$total*100,"\n";
#print "Merged_seq_with_high_quality: $flag\t",$flag/$total*100,"\n";
#print "Unmerged_seq: $unmerge_n\t",$unmerge_n/$total*100,"\n";
#print "Merged_but_filter_by_low_quality: $filter_n\t",$filter_n/$total*100,"\n";
print "Unique_seq_num: ",scalar keys %combine_seq,"\n";

