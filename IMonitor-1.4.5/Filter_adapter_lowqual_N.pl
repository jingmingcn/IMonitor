#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

=head1	Name
	Filter_adapter_lowqual_N.pl	filter adapter pollution, low quality base/sequence, too much "N" base
	
=head1	Usage
	perl Filter_low_N_Qual.pl <*_1.fq(.gz)> <*_2.fq(.gz)> <read_length> <out_dir> <output_file_index> <seqtype> [<adapter1.list> <adapter2.list>] [options]

=head1	Options:

	-F      <Flag> filter the adapter pollution
	-N	<FLOAT> filter the read with N content [0.05]
	-1	<INT> minimum length of read1 retained [read length/2+10]
	-2	<INT> minimum length of read2 retained [read length/2]
	-Q	<INT> filter the low quality base of read's end [10]
	-L	<INT> adapter pollution in last Xbp of read, then cut the adapter sequence [50]
	-v	<INT> used to calculate the base quality [64]

=head1	Instruction:
	Conduct the sequences with the conditions as follow:
	1. filter read with some "N" bases
	2. filter the adapter pollution: if the adapter in the last Xbp, cut the adapter,otherwise, filter the sequence
	3. cut read's end base with low quality
	4. restrict the retained read length


=cut

my ($N_f_r,$read1_min_l , $read2_min_l , $qual_cut, $FLAG, $last_len, $cal_qual_v);
$FLAG = 0;

GetOptions(
		"F" => \$FLAG,
		"N:i" => \$N_f_r,
		"1:i" => \$read1_min_l,
		"2:i" => \$read2_min_l,
		"Q:i"=> \$qual_cut,
		"L:i", => \$last_len,
		"v:i" => \$cal_qual_v
	  );
die `pod2text $0` unless(@ARGV == 6 || @ARGV == 8);
my ($in_file1 , $in_file2 , $read_length, $dir , $f_index, $seqtype, $adap_1 , $adap_2) = @ARGV;
$dir =~s/\/$//;


# default values for parameters
$N_f_r||=0.05;
$read1_min_l||=$read_length/2+10;
$read2_min_l||=$read_length/2;
$qual_cut||=10;
$last_len||=$read_length/2;
$cal_qual_v=64 unless(defined($cal_qual_v));


# 1. read the adapter list files and instore in two hash
my (%adapter_list1,%adapter_list2);
if($FLAG){
	die "Need to filter adapter pollution, but the file\"adapter.list.gz\" were not found!\n" unless(defined($adap_1) && defined($adap_1));
	if($adap_1=~/\.gz/){open I, "gzip -dc $adap_1|" or die;}
	else{open I, "$adap_1" or die;}
	<I>;
	while(<I>){
		chomp;
		my ($read_name,$name2,$left,$right);
		if($seqtype){
			($read_name,$name2,$left,$right)=(split /\s+/,$_)[0,1,3,4];
			$read_name = "$read_name-$name2";
		}else{
			($read_name,$left,$right)=(split /\s+/,$_)[0,2,3];
		}
		$adapter_list1{$read_name}=[($left,$right)];
	}
	close I;
	if($adap_2=~/\.gz/){open I, "gzip -dc $adap_2|" or die;}
	else{open I, "$adap_2" or die;}
	<I>;
	while(<I>){
		chomp;
		my ($read_name,$name2,$left,$right);
		if($seqtype){
			($read_name,$name2,$left,$right)=(split /\s+/,$_)[0,1,3,4];
			$read_name = "$read_name-$name2";
		}else{
			($read_name,$left,$right)=(split /\s+/,$_)[0,2,3];
		}
		$adapter_list2{$read_name}=[($left,$right)];
	}
	close I;
}

# 2. read the fq.gz files
if($in_file1=~/\.gz/){open IN1 , "gzip -dc $in_file1|" or die;}
else{open IN1 , "$in_file1" or die;}
if($in_file2=~/\.gz/){open IN2 , "gzip -dc $in_file2|" or die;}
else{open IN2 , "$in_file2" or die;}

open O1 , "|gzip>$dir/${f_index}_filter_1.fq.gz" or die;
open O2 , "|gzip>$dir/${f_index}_filter_2.fq.gz" or die;

my $read_sum = 0;
my $retain_num = 0;
my $retain_read1_l_ave = 0;
my $retain_read2_l_ave = 0;
my $read1_cut_n = 0;
my $read2_cut_n = 0;
my $filter_N_n = 0;
my $filter_Q_n = 0;
my $adapter_reads = 0;

while(<IN1>)
{
	chomp;
	my $id_1 = $_;
	$id_1 =~s/\s+/-/;
	chomp(my $seq_1 = <IN1>);
	chomp(my $flag_1 = <IN1>);
	chomp(my $qual_1 = <IN1>);

	chomp(my $id_2 = <IN2>);
	$id_2 =~s/\s+/-/;
	chomp(my $seq_2 = <IN2>);
	chomp(my $flag_2 = <IN2>);
	chomp(my $qual_2 = <IN2>);
	$read_sum++;

        #       0. filter reads with "N" number
	my $N_1_n = $seq_1 =~ tr/Nn/Nn/;
	my $N_2_n = $seq_2 =~ tr/Nn/Nn/;
	if($N_1_n >= $N_f_r*$read_length || $N_2_n >= $N_f_r*$read_length){
		$filter_N_n++;
		next;
	}



	#	1. filter the adapter pollution
	if($FLAG)
	{
		my $read_id1=$id_1;$read_id1=~s/^@//;
		my $read_id2=$id_2;$read_id2=~s/^@//;
		if(exists $adapter_list1{$read_id1})
		{
			if($adapter_list1{$read_id1}->[0]<=($read_length-$last_len)){
				$adapter_reads++;
				next;
			}
			else{
				$seq_1=substr($seq_1,0,$adapter_list1{$read_id1}->[0]);
				$qual_1=substr($qual_1,0,$adapter_list1{$read_id1}->[0]);
			}
		}

		if(exists $adapter_list2{$read_id2})
		{
			if($adapter_list2{$read_id2}->[0]<=($read_length-$last_len)){
				$adapter_reads++;
				next;
			}
			else{
				$seq_2=substr($seq_2,0,$adapter_list2{$read_id2}->[0]);
				$qual_2=substr($qual_2,0,$adapter_list2{$read_id2}->[0]);
			}
		}

	}


	my ($r_l_1 ,$r_l_2)= (length($seq_1) , length($seq_2));
	

	#	2. filter read1 or cut end bases of read1
	my @seq_1_b = split // , $seq_1;
	my @qual_1_b = split // , $qual_1;

	
	for(my $i=$#seq_1_b ; $i>=$read1_min_l-1 ; $i--){
		# 2-1. cut the end base with quality lower than $qual_cut
		if(ord($qual_1_b[$i])-$cal_qual_v <= $qual_cut){
			pop @seq_1_b;
			pop @qual_1_b;
			next;
		}
		else{
			last;
		}
	}
	# 2-3. retained sequence length must be >= $read1_min_l
	$seq_1 = join "" ,@seq_1_b;
	$qual_1 = join "" ,@qual_1_b;
	next if(length($seq_1)<$read1_min_l);

	#	3. filter read2 or cut end bases of read2
	my @seq_2_b = split // , $seq_2;
	my @qual_2_b = split // , $qual_2;

	for(my $i=$#seq_2_b ; $i>=$read2_min_l-1 ; $i--){
		# 2-1. cut the end base with quality lower than $qual_cut
		if(ord($qual_2_b[$i])-$cal_qual_v <= $qual_cut){
			pop @seq_2_b;
			pop @qual_2_b;
			next;
		}
		else{
			last;
		}
	}
	# 2-3. retained sequence length must be >= $read2_min_l
	$seq_2 = join "" ,@seq_2_b;
	$qual_2 = join "" ,@qual_2_b;
	next if(length($seq_2)<$read2_min_l);

	#	print out
	$retain_num++;
	$retain_read1_l_ave = $retain_read1_l_ave+length($seq_1);
	$retain_read2_l_ave = $retain_read2_l_ave+length($seq_2);
	$read1_cut_n++ if(length($seq_1)<$r_l_1);
	$read2_cut_n++ if(length($seq_2)<$r_l_2);

	print O1 "$id_1\n$seq_1\n$flag_1\n$qual_1\n";
	print O2 "$id_2\n$seq_2\n$flag_2\n$qual_2\n";
	
}
close IN1;
close IN2;



# print out statistic information
$filter_Q_n = $read_sum - $retain_num - $filter_N_n - $adapter_reads;
($retain_read1_l_ave , $retain_read2_l_ave) = ($retain_read1_l_ave/$retain_num , $retain_read2_l_ave/$retain_num);
print "All_read_number(PE=1):	$read_sum\n";
print "retained_read_num: $retain_num\t",$retain_num/$read_sum*100,"\n";
print "filter_with_N_num: $filter_N_n\n";
print "filter_with_Adapter: $adapter_reads\n";
print "filter_with_Length: $filter_Q_n\n";
print "read1's_end_with_triming: $read1_cut_n\n";
print "read2's_end_with_triming: $read2_cut_n\n";
print "average_length_of_retain_read1: $retain_read1_l_ave\n";
print "average_length_of_retain_read2: $retain_read2_l_ave\n";


