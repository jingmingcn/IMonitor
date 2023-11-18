#!/usr/bin/perl -w
use strict;

=head1 Usage
	perl Merge_filter_fqtofa.pl  <cope.merge.fq> <java.merge.fq.gz> <fqMerging.log> <filter-Q> <filter-rate> <out.fa.gz> <change_id.back.gz> <insert-size.distribution.gz><quality_system>
	<filter-rate: the rate of low quality bases to filter>
	<out.fa.gz>/<change_id.back.gz>: must be the suffix ".gz"

=head1 Description
	1. filter the mereged sequence with too much low quality bases
	2. change fq to fa
	3. insert size length distribution statistic

=cut

die `pod2text $0` unless(@ARGV==9);


my ($cope_f, $java_f, $log, $Q, $rate , $out , $out_id , $out_insert,$cal_qual_v) = @ARGV;

open I, "$cope_f" or die;

open OUT, "|gzip >$out" or die;
open OD, "|gzip >$out_id" or die;
#open OF, "|gzip >$out_filter" or die;

my %insertsize_len;
my $java_meger_n = 0;
my $filter_n = 0;
my $flag = 0;
my $sum = 0;
my $flag_c = 0;
my %combine_seq;

# 1. read the cope-mereged file
COPE:while(<I>)
{
	chomp;
	my $id = $_;
	chomp(my $seq = <I>);
	chomp(my $f = <I>);
	chomp(my $qual = <I>);

	# 1.1 insert size length 
	my $len = length($seq);
	$insertsize_len{$len}++;
	$sum++;

	# 1.2 filter the low quality sequence	
	my $low_n_a = $len*$rate;
	my @all_q = split // , $qual;
	my $low_n = 0;
	for(@all_q){
		if($low_n > $low_n_a){
			$filter_n++;
#			print OF "$id\n$seq\n$f\n$qual\n";
			next COPE;
		}
		$low_n++ if(ord($_)-$cal_qual_v <= $Q);
	}

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

open I, "gzip -dc $java_f|" or die;

# 2. read the jave-mereged file
JAVA:while(<I>)
{
	$java_meger_n++;
	chomp;
        my $id = $_;
        chomp(my $seq = <I>);
        chomp(my $f = <I>);
        chomp(my $qual = <I>);

        # 1.1 insert size length
        my $len = length($seq);
        $insertsize_len{$len}++;
	$sum++;

        # 1.2 filter the low quality sequence
        my $low_n_a = $len*$rate;
        my @all_q = split // , $qual;
        my $low_n = 0;
        for(@all_q){
                if($low_n > $low_n_a){
                        $filter_n++;
#                        print OF "$id\n$seq\n$f\n$qual\n";
                        next JAVA;
                }
                $low_n++ if(ord($_)-$cal_qual_v <= $Q);
        }

        # 1.3 change fq to fa
        $flag++;
        my $new_id = "cp$flag";
        $id = (split)[0];
	my $com_id;
	if(exists $combine_seq{$seq}){
		my ($d,$n) = split /:/,$combine_seq{$seq};
		$n++;
		$combine_seq{$seq} = "$d:$n";
		$com_id = $d;
	}else{
		$flag_c++;
		$com_id = "cp${flag_c}c";
		$combine_seq{$seq} = "$com_id:1";
	}

        print OD "$id\t$new_id\t$len\t$com_id\n";
#        print OUT ">$new_id\n$seq\n";
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

my $tmp = `awk 'NR==3{print \$1"\t"\$2"\t"\$3}' $log`;
chomp($tmp);
my ($total, $cope_n , $cope_r) = split /\t/,$tmp;
my $unmerge_n = $total-$cope_n-$java_meger_n;

print "Total_seq: $total\n";
print "Merged_seq_raw: ",$cope_n+$java_meger_n,"\t",($cope_n+$java_meger_n)/$total*100,"\n";
print "\tCOPE_merged_seq_raw: $cope_n\t",$cope_n/$total*100,"\n";
print "\tJava_merged_seq_raw: $java_meger_n\t",$java_meger_n/$total*100,"\n";
print "Merged_seq_with_high_quality: $flag\t",$flag/$total*100,"\n";
print "Unmerged_seq: $unmerge_n\t",$unmerge_n/$total*100,"\n";
print "Merged_but_filter_by_low_quality: $filter_n\t",$filter_n/$total*100,"\n";
print "Unique_seq_num: ",scalar keys %combine_seq,"\n";

