#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin $Script);

=head1 Name 
	Correct_seq_pcr_error.pl 

=head1 Description
	

=head1 Usage
	
	perl Correct_seq_pcr_error.pl -i -d -o [options]

	-i	<FILE>	input file
	-Q 	<INT>	low quality cutoff[20]
	-M1	<INT>	low quality bases number allowed [5]
	-M2	<INT>	PCR orror number[3]
	-R	<INT>	the multiple of abundance allowed when consider the PCR error[5]
	-d	<FILE>	filter sequence (*.gz) 
	-o	<FILE>	output file (*.gz)
	-c		logical value, correct only for CDR3
	-n	<INT>   sequence abundance, if more than this, sequencing error will not be considered[50]
	-v      <I> used to calculate the base quality [64]


=cut

my ($low_q ,$low_q_n , $abundance_multi , $pcr_e_n , $in , $left_f , $out, $cdr3_flag, $maximum_abund , $q_system);
GetOptions(
		"i=s"=>\$in,
		"Q:i"=>\$low_q,
		"M1:i"=>\$low_q_n,
		"R:i"=>\$abundance_multi,
		"M2:i"=>\$pcr_e_n,
		"d=s"=>\$left_f,
		"o=s"=>\$out,
		"c"=>\$cdr3_flag,
		"n:i"=>\$maximum_abund,
		"v:i"=>\$q_system
	  );

die `pod2text $0` if(!$in && !$left_f && !$out);

# default paremeters
$low_q=20 unless(defined($low_q));
$low_q_n=5 unless(defined($low_q_n));
$abundance_multi=5 unless(defined($abundance_multi));
$pcr_e_n=3 unless(defined($pcr_e_n));
$maximum_abund= 50 unless(defined($maximum_abund));
$q_system=64 unless(defined($q_system));

if($in=~/\.gz$/){
	open I, "gzip -dc $in|" or die;
}else{
	open I, "$in" or die;
}
open D, "|gzip >$left_f" or die;
open O, "|gzip >$out" or die;



my %substitute;
$substitute{"A"} = [("C","G","T")];
$substitute{"C"} = [("A","G","T")];
$substitute{"G"} = [("A","C","T")];
$substitute{"T"} = [("A","G","C")];
$substitute{"N"} = [("A","G","C","T")];


my $solexa_q;
open A, "$Bin/ASCII.txt" or die;
while(<A>)
{
	chomp;
	my @line = split;
	if($line[2]-$q_system >=0 && $line[2]-$q_system  < $low_q){
		if(defined($solexa_q)){
			$solexa_q .= $line[4];
		}else{
			$solexa_q = $line[4];
		}
	}
}
close A;

#print "$solexa_q\n";

#------------------	1. read input file && filter the sequence more than $low_q_n low quality 
my %core_seq;
my %low_qual_seq;
my %core_seq_backup;
my %fir_core_record;
my %low_qula_abund;

my ($all_n , $first_core , $low_q_f , $low_q_c , $low_q_unmap , $second_core , $second_c) = (0,0,0,0,0,0,0);

#my $t =localtime(time);
#print "read file:$t\n";


while(<I>)
{
	$all_n++;
	chomp;
	my $raw = $_;
	my @line = split;
	my ($l_n , $lowq_s)= &find_low_qual($line[33],$line[32]); # low quality base count
	# 1) too much low qualities
	if($l_n > $low_q_n){
		$low_q_f++;
		print D "$raw\n";
	}
	# 2) core sequence
	elsif($l_n == 0){
		$first_core++;
		my $v_subf = (split /[*.]/,$line[0])[0];
		my $j_subf = (split /[*.]/,$line[20])[0];
		$core_seq{"$v_subf$j_subf"}{$line[32]}++;
		$core_seq_backup{$line[30]} = 1;
		$fir_core_record{"$v_subf$j_subf\t$line[32]"} = $raw;
	}
	# 3) low quality <= $low_q+n
	else{
		$low_qual_seq{$raw} = $lowq_s;
		my $v_subf = (split /[*.]/,$line[0])[0];
		my $j_subf = (split /[*.]/,$line[20])[0];
		$low_qula_abund{"$v_subf$j_subf\t$line[32]"}++;
	}

}
close I;

#$t =localtime(time);
#print "finished read:$t\n";

#---------------------	2. map the low quality sequenct to core sequence	----------------#
my %correct_low_seq;

for my $r(keys %low_qual_seq)
{
	my @line = split /\s+/,$r;
	my $v_subf = (split /[*.]/,$line[0])[0];
	my $j_subf = (split /[*.]/,$line[20])[0];

	# 1) sequence is the same with one of core sequence
	if(exists $core_seq{"$v_subf$j_subf"}{$line[32]} && $low_qula_abund{"$v_subf$j_subf\t$line[32]"}>=$maximum_abund)
	{
		$low_q_c++;
		$core_seq{"$v_subf$j_subf"}{$line[32]}++;
		$core_seq_backup{$line[30]} = 1;
		next;
	}
	
	# 2) map the sequence which mismatch allowed mapp to core sequence
	my ($flag,$new) = &construct_matrix($low_qual_seq{$r},$line[32],"$v_subf$j_subf");# construct a matrix sequence with mismatch
	if($flag){
		$low_q_c++;
		$core_seq{"$v_subf$j_subf"}{$new}++;
		$core_seq_backup{$line[30]} = 1;
		$correct_low_seq{$line[30]} = "$v_subf$j_subf\t$new";
	}
	else{
		$low_q_unmap++;
		print D "$r\n";
	}
}
%low_qual_seq = ();
%low_qula_abund = ();
$second_core = $first_core + $low_q_c;

# $t =localtime(time);
#print "finish Seq error:$t\n";


#----------------- 3. remove the PCR error according sequence's abundance  --------------------------#
my %correct_pcr_seq;

#my %sub_max;
#for my $s1 (keys %core_seq){
#	my $t = (sort {$b<=>$a} values %{$core_seq{$s1}})[0];
#	$sub_max{$s1} = $t;
#}


#for my $subfam1 (keys %core_seq)
#{
#	for my $seq1 (sort {$core_seq{$subfam1}{$a}<=>$core_seq{$subfam1}{$b}} keys %{$core_seq{$subfam1}})
#	{
#	last if($core_seq{$subfam1}{$seq1}*$abundance_multi > $sub_max{$subfam1} );

#	my ($min_mis , $max_abund , $corr_seq) = (100 , 0);
#	my $flag = 0;
#	for my $seq2 (sort {$core_seq{$subfam1}{$b}<=>$core_seq{$subfam1}{$a}} keys %{$core_seq{$subfam1}})
#	{
#		last if($core_seq{$subfam1}{$seq2}<= $core_seq{$subfam1}{$seq1}*$abundance_multi); # if the abundane less than $abundance_multi, need not to consider

#		next if(length($seq1) ne length($seq2));# only consider the sequences with same length
		
#		my $mis_n = &seq_align($seq1 , $seq2);# calculate the mismatch of two sequences
#		if($mis_n <= $pcr_e_n)# probale the PCR error
#		{
#			$flag = 1;
#			if($mis_n < $min_mis){
#				$min_mis = $mis_n;
#				$max_abund = $core_seq{$subfam1}{$seq2};
#				$corr_seq = "$subfam1\t$seq2";
#			}
#			elsif($mis_n == $min_mis){
#				$max_abund = $core_seq{$subfam1}{$seq2} if($max_abund<$core_seq{$subfam1}{$seq2});
#				$corr_seq = "$subfam1\t$seq2";
#			}
#		}

#	}
#	if($flag){
#		$correct_pcr_seq{"$subfam1\t$seq1"} = $corr_seq;
#		$second_c += $core_seq{$subfam1}{$seq1};
#	}
#	}
#}

#$t =localtime(time);
#print "finish pcr err:$t\n";

#-----------------------------	4. print out the effective data	----------------------------------------#
if($in=~/\.gz$/){
	open I, "gzip -dc $in|" or die;
}else{
	open I, "$in" or die;
}
while(<I>)
{
	chomp;
	my @line = split;
	next if(!exists $core_seq_backup{$line[30]});
	
	my $v_subf = (split /[*.]/,$line[0])[0];
	my $j_subf = (split /[*.]/,$line[20])[0];
	if(exists $correct_low_seq{$line[30]})
	{
		my $new_r;
		if(exists $correct_pcr_seq{"$v_subf$j_subf\t$line[32]"}){
			$new_r = $fir_core_record{$correct_pcr_seq{"$v_subf$j_subf\t$line[32]"}};
		}
		else{
			$new_r = $fir_core_record{$correct_low_seq{$line[30]}};
		}

		if($cdr3_flag)# only CDR3 correction
		{
#			my $cdr3_ref_start = (split /\./,$line[0])[4];
#			my $cdr3_start = $line[5]-($line[7]-$cdr3_ref_start);
			my $new_cdr3 = (split /\s+/,$new_r)[32];
#			substr($line[31],$cdr3_start-1,length($new_cdr3)) = $new_cdr3;
			substr($line[31],index($line[31],$line[32]),length($new_cdr3)) = $new_cdr3;
			$new_r = join "\t", @line;

		}

		my @l_new_r = split /\s+/,$new_r;
		pop @l_new_r;pop @l_new_r;pop @l_new_r;
		$l_new_r[30] = $line[30];
		print O join "\t" , @l_new_r;print O "\n";
	}
	elsif(exists $correct_pcr_seq{"$v_subf$j_subf\t$line[32]"})
	{
		my $new_r = $fir_core_record{$correct_pcr_seq{"$v_subf$j_subf\t$line[32]"}};
		if($cdr3_flag)# only CDR3 correction
		{
#			my $cdr3_ref_start = (split /\./,$line[0])[4];
#			my $cdr3_start = $line[5]-($line[7]-$cdr3_ref_start);
			my $new_cdr3 = (split /\s+/,$new_r)[32];
#			substr($line[31],$cdr3_start-1,length($new_cdr3)) = $new_cdr3;
			substr($line[31],index($line[31],$line[32]),length($new_cdr3)) = $new_cdr3;
			$new_r = join "\t", @line;
		}
		

		my @l_new_r = split /\s+/,$new_r;
		pop @l_new_r;pop @l_new_r;pop @l_new_r;
		$l_new_r[30] = $line[30];
		print O join "\t" , @l_new_r;print O "\n";
	}
	else
	{
		pop @line;pop @line;pop @line;
		print O join "\t",@line;print O "\n";
	}
}
close I;
close O;
print "all_raw_seq: $all_n\n";
print "high_quality_seq: $first_core\n";
print "low_quality_correct: $low_q_c\n";
print "low_quality_filter: $low_q_f\n";
print "low_quality_unmerge: $low_q_unmap\n\n";
print "new_core_sequence: $second_core\n";
print "PCR_sequence_correct: $second_c\n";


#--- calculate the number of low quality base  ---#
sub find_low_qual
{
	my ($q,$s) = @_;
	my $l_n = 0;
#	$l_n = $q=~tr/@ABCDEFGHIJKLMNOPQRS/@ABCDEFGHIJKLMNOPQRS/;
	eval "\$l_n = \$q=~tr/$solexa_q/$solexa_q/";
	return $l_n if($l_n > $low_q_n);

	my @q_all = split //, $q;
	my @s_all = split //,$s;
	my $low_q_all;
	for(my $i=0; $i<=$#q_all ; $i++){
		if(ord($q_all[$i])-$q_system < $low_q){
			$l_n++;
			if(!defined($low_q_all)){
				$low_q_all = "$i:$s_all[$i]";
			}
			else{
				$low_q_all .= ":$i:$s_all[$i]";
			}
		}
		last if($l_n > $low_q_n);
	}
	if($l_n > $low_q_n){
		return $l_n;
	}else{
		return ($l_n , $low_q_all);
	}

}

#---	construct a matrix sequence with mismatch
# store in a array 
# less mismatch 's sequence in the former of array and more mismatch's sequences in the back
#


sub construct_matrix
{
	my ($low_f , $r_s, $vj) = @_;
	my @matrix_m;
	my %low_f_a = split /:/,$low_f;
	push @matrix_m , $r_s;

	for(my $i=0 ; $i<scalar keys %low_f_a ; $i++)
	{
		for my $j(keys %low_f_a){
			my @bases;
			my $b = $low_f_a{$j};
			@bases = @{$substitute{$b}};

			for my $old_s (@matrix_m)
			{
				for(@bases){
					my $new = $old_s; 
					substr($new , $j , 1)= $_;
					if(exists $core_seq{"$vj"}{$new}){
						return (1,$new);
					}
				}
			}
		}
	}
	return (0,$r_s);
}

#---------- seq alignment	---------------------#
sub seq_align
{
	my ($s1 , $s2) = @_;
	my $sub_len = int(length($s1)/($pcr_e_n+1));
	my ($s1_n , $s2_n);
	my $mis_n = 0;
	my $flag = 0;
	for(my $i=0 ; $i<$pcr_e_n ; $i++)
	{
		my $s1_sub = substr($s1, $i*$sub_len , $sub_len);
		my $s2_sub = substr($s2, $i*$sub_len , $sub_len);
		if($s1_sub eq $s2_sub){
			$flag = 1;
		}
		else{
			$s1_n .= $s1_sub;
			$s2_n .= $s2_sub;
		}
	}
	my $s1_sub = substr($s1, $pcr_e_n*$sub_len);
	my $s2_sub = substr($s2, $pcr_e_n*$sub_len);
	if($s1_sub eq $s2_sub){
		$flag = 1;
	}
	else{
		$s1_n .= $s1_sub;
		$s2_n .= $s2_sub;
	}

	if($flag){
		my @s1_a = split //,$s1_n;
		my @s2_a = split //,$s2_n;
		for(my $j=0 ; $j<=$#s1_a ; $j++){
			$mis_n++ if($s1_a[$j] ne $s2_a[$j]);
			last if($mis_n > $pcr_e_n);
		}
		return $mis_n;
	}
	else{
		return ($pcr_e_n+1);		
	}
}

