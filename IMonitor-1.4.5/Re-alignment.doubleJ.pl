#!/usr/bin/perl -w
use Getopt::Long;
use Data::Dumper;
use strict;

#	Author:	 zhangwei3@genomics.cn
# Modify 2012/12/21:
#       conduct the VJ overlap region: cut the J align part which located in the V aligned region
#       old: conduct everty aligned record's edge , then select the best align ==> new: select the best align ,then conduct the best align record's edge
#
# Modify 2012/12/24:
#	change the method of caculating score and identity,
#		old: all ref sequence length, new: the CDR3 part's edge ensured by ref seq, the other part's edge ensured by sequencing seq
# Modify 2013/05/07:
#	chang the re-align method:
#	1. caculate score and identity: non-CDR3 region use global alignment, CDR3 region use local alignent. 
#	   1) non-CDR3 region : extended to the whole sequence according to BLAST result(V: 5-> ; J: 3->)
#	   2) CDR3 region: extend with restricted mismatch, and remove the all mismatch of end(e.g, we allow 2 mismatch,find the 3th mismatch and stop 3th-mismathc-1 postion, then re-back and remove the mismatch of end sequence)
#	   3) according to 1) and 2) 's alignment part, caculate identity and sorce
#	2. selcet the highest score as the best hit. if there are several highest score records, then select the least deletion number 's record(VJ). for D ,then select the most high identity
#
# Modify 2013/05/11:
#	1. add the VDJ alignemt statistic and screen output
# 	2. use -iv -ij paremeters instead -I
# Modify 2015/08:
#	the part(seq) exceeded 	ref end do not consier for align-length,mismatch,score,identity.
# Modify 2019/3:
# 	my $cdr3_p = (split /\./,(split /\s+/ , $record)[1])[4];  change to my $cdr3_p = (split /\./,(split /\s+/ , $record)[1])[-2];
#	Lines number: 472,712

my $usage=<<USAGE;

        Usage:
        perl $0 <fa.gz><blast.V.m9.gz><blast.J.m9.gz><out_dir_file><sample_name><VDJ.converted.fa><type>

Options:
	-F	<file>	the file of D blast
        -v	<INT>	the mismatch of V during CDR3 region [TRA/TRB:0,IGH:2,IGKL:7]
	-d	<INT>	the mismatch of D during CDR3 region [TRB:0,IGH:4]
	-j 	<INT>	the mismatch of J during CDR3 region [TRB:0,IGH:2]
        -vif	<FLOAT>	V align identity filter [80]
	-jif	<FLOAT>	J align identity filter [80]
	-r	<INT>	VJ align minmun length restrict [6]
	-l	<INT>	the length restrict of D align [4]
	-m	<INT>	match penalty [5]
	-p	<INT>	mismtch penalty [-4]

USAGE

my ($fa , $v_f , $j_f , $out_dir , $sample , $ref_f , $type) = @ARGV;
my ($d_f,$v_mismatch,$d_mismatch,$j_mismatch,$v_identity_filter,$j_identity_filter,$D_len_restrict,$match_s,$mismatch_s,$VJ_len_restrict);

GetOptions(
		'F:s'=>\$d_f,
		'v:i'=>\$v_mismatch,
		'd:i'=>\$d_mismatch,
		'j:i'=>\$j_mismatch,
		'vif:f'=>\$v_identity_filter,
		'jif:f'=>\$j_identity_filter,
		'r:i'=>\$VJ_len_restrict,
		'l:s'=>\$D_len_restrict,
		'm:i'=>\$match_s,
		'p:i'=>\$mismatch_s
		);

die $usage unless($fa && $v_f && $j_f && $out_dir && $sample && $ref_f && $type); 

$v_identity_filter =80 unless(defined($v_identity_filter));
$j_identity_filter=80 unless(defined($j_identity_filter));
$D_len_restrict=4 unless(defined($D_len_restrict));
$match_s=5 unless(defined($match_s));
$mismatch_s=-4 unless(defined($mismatch_s));
$VJ_len_restrict=6 unless(defined($VJ_len_restrict));

my $D_flag = 0;
$D_flag = 1 if(defined($d_f));# need to analyze D gene

# initialization $v_mismatch and $j_mismatch
unless(defined($v_mismatch)){
	$v_mismatch = 0 if($type eq "TRB" || $type eq "TRA" || $type eq "TRD" || $type eq "TRG" || $type eq "TRAB" || $type eq "TRDG" || $type eq "TR");
	$v_mismatch = 2 if($type eq "IGH" || $type eq "IG");
	$v_mismatch = 7 if($type eq "IGKL");
}
unless(defined($j_mismatch)){
	$j_mismatch = 0 if($type eq "TRB" || $type eq "TRA" || $type eq "TRD" || $type eq "TRG" || $type eq "TRAB" || $type eq "TRDG" || $type eq "TR");
	$j_mismatch = 2 if($type eq "IGH" || $type eq "IG");
	$j_mismatch = 7 if($type eq "IGKL");
}
unless(defined($d_mismatch)){
	$d_mismatch = 0 if($type eq "TRB" || $type eq "TRD" || $type eq "TRAB" || $type eq "TRDG" || $type eq "TR");
	$d_mismatch = 4 if($type eq "IGH" || $type eq "IG");
}

open OUT , ">$out_dir" or die;

#--------------------------	read ref file	------------------
my %ref;# key: seq-id , value: seq

open I, "$ref_f" or die;
$/ = ">";<I>;
while(<I>)
{
	chomp;
	my @line=split /\n+/,$_;
	my $tile=shift @line;
	my $seq="@line";
	$seq=~s/\s+//g;
	$ref{$tile} = $seq;
}
$/ = "\n";
close I;

#------------------------- read merge-fa file	--------------------------------#

my ($pre_v , $pre_j , $pre_d);

if($fa=~/\.gz/){open I, "gzip -dc $fa|" or die;}
else{open I,"$fa" or die;}

if($v_f=~/\.gz/){open V, "gzip -dc $v_f|" or die;}
else{open V, "$v_f" or die;}
chomp($pre_v = <V>);


if($j_f=~/\.gz/){open J, "gzip -dc $j_f|" or die;}
else{open J,"$j_f" or die;}
chomp($pre_j=<J>);

if($D_flag)
{
	if($d_f=~/\.gz/){open D, "gzip -dc $d_f|" or die;}
	else{open D, "$d_f" or die;}
	chomp($pre_d = <D>);
}
# read mergeed_fa file	
my ($fa_num , $V_align_n , $D_align_n , $J_align_n) = (0,0,0,0);
my %vdj_n_stat;

while(<I>)
{
	chomp;
	my @L_a = split;
	my $seq_id = $L_a[30];#print "$seq_id\n";
	#	$seq_id =~ s/^>//;
	#	chomp(my $seq = <I>);
	my $seq = $L_a[31];
	my $abund_n = (split /:/,$seq_id)[1];
	$fa_num+=$abund_n;
	
	my $double_J = 0;
	my ($final_V,$final_D,$final_J);
	# judge whether it need to realign again
	
	$final_V = join "\t", @L_a[0..9];
	$final_D = join "\t", @L_a[10..19];
	$final_J = join "\t", @L_a[20..29];
	if($L_a[0] ne "NA" && $L_a[20] ne "NA")
	{
		my @N = sort {$a<=>$b}($L_a[4],$L_a[5],$L_a[24],$L_a[25]);
		my $interval = $N[2]-$N[1];
		if($interval>100){
			$double_J = 1;
		}
	}
	

	my ($v_end_pos , $j_start_pos , $strand_v);
	my (%v_align_record , %d_align_record , %j_align_record); # store one sequence 's align records
	

	# -----------------------	conduct V gene	---------------------

	if($pre_v eq "EOF" || (split /\s+/, $pre_v)[0] ne $seq_id) # the sequence without align record
	{
		print OUT "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
		$v_end_pos = "NA";
		$strand_v = "NA";
	}
	else
	{
		# read V align file
	 	$v_align_record{$pre_v} = 1;
		
		while(<V>)
		{
			chomp;
			if((split /\s+/, $_)[0] eq $seq_id)# store records
			{
				$v_align_record{$_} = 1;
			}
			else # read the next seq record
			{
				my $best_v ;
				if($double_J)# re-alignment
				{
					($best_v , $v_end_pos , $strand_v)= &best_align_vj("V" , $seq , \%v_align_record); # select the best V align
					$best_v = join "\t" , (split /\s+/,$best_v)[1,2,3,4,6,7,8,9,10,11];
				}else{
					$best_v = $final_V;
				}
				print OUT "$best_v"; # print out
				if((split /\t/ , $best_v)[0] ne "NA"){
					$vdj_n_stat{$seq_id}=1;
					$V_align_n+=$abund_n;
				}
				$pre_v = $_;
				%v_align_record = ();
				last;
			}
		}
		
		if((scalar keys %v_align_record) != 0) # the last records of V align file
		{
			my $best_v ;
			if($double_J)# re-alignment
			{
				($best_v , $v_end_pos , $strand_v)= &best_align_vj("V" , $seq , \%v_align_record); # select the best V align
				$best_v = join "\t" , (split /\s+/,$best_v)[1,2,3,4,6,7,8,9,10,11];
			}else{
				$best_v = $final_V;
			}
			print OUT "$best_v"; # print out
			if((split /\t/ ,$best_v)[0] ne "NA"){
				$vdj_n_stat{$seq_id}=1;
				$V_align_n+=$abund_n;
			}
			$pre_v = "EOF";
			%v_align_record = ();
		}
	}


	# -----------------------       conduct J gene  ---------------------

	my $best_j ;
	if($pre_j eq "EOF" || (split /\s+/, $pre_j)[0] ne $seq_id) # the sequence without align record
	{
		$best_j = "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
		$j_start_pos = "NA";
	}
	else
	{
		# read V align file
	 	$j_align_record{$pre_j} = 1;
		while(<J>)
		{
			chomp;
			if((split /\s+/, $_)[0] eq $seq_id)# store records
			{
				$j_align_record{$_} = 1;
			}
			else# read the next seq record
			{
				if($double_J)# re-alignment
				{
					&Modify_J_records(\%j_align_record,$strand_v);
					($best_j , $j_start_pos)= &best_align_vj("J" , $seq , \%j_align_record , $v_end_pos , $strand_v); # select the best V align 
					$best_j= join "\t" , (split /\s+/,$best_j)[1,2,3,4,6,7,8,9,10,11];
				}else{
					$best_j = $final_J;
				}
				if((split /\t/ ,$best_j)[0] ne "NA"){
					$vdj_n_stat{$seq_id}++ if(exists $vdj_n_stat{$seq_id});
					$J_align_n+=$abund_n;
				}
				$pre_j = $_;
				%j_align_record = ();
				last;
			}
		}
		
		if((scalar keys %j_align_record) != 0) # the last records of V align file
		{
			if($double_J)# re-alignment
			{
				&Modify_J_records(\%j_align_record,$strand_v);
				($best_j , $j_start_pos)= &best_align_vj("J" , $seq , \%j_align_record , $v_end_pos , $strand_v); # select the best V align
				$best_j= join "\t" , (split /\s+/,$best_j)[1,2,3,4,6,7,8,9,10,11];
			}else{
				$best_j = $final_J;
			}
			if((split /\t/, $best_j)[0] ne "NA"){
				$vdj_n_stat{$seq_id}++ if(exists $vdj_n_stat{$seq_id});
				$J_align_n+=$abund_n;
			}
			%j_align_record = ();
			$pre_j = "EOF";
		}
	}


	# -----------------------       conduct D gene  ---------------------
	unless($D_flag) # no D gene 
	{
		print OUT "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$best_j";
		print OUT "\t$seq_id\t$seq\n";
		next;
	}
	if($pre_d eq "EOF" || (split /\s+/, $pre_d)[0] ne $seq_id) # the sequence without align record
	{
		print OUT "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$best_j";
	}
	else
	{
		# read D align file
	 	$d_align_record{$pre_d} = 1;	
		while(<D>)
		{
			chomp;
			if((split /\s+/, $_)[0] eq $seq_id)# store records
			{
				$d_align_record{$_} = 1;
			}
			else # read the next align record
			{
				my $best_d ;
				if($double_J)# re-alignment
				{
					$best_d = &best_align_d($seq , \%d_align_record , $v_end_pos , $j_start_pos , $strand_v); # select the best D align 
					$best_d= join "\t" , (split /\s+/,$best_d)[1,2,3,4,6,7,8,9,10,11];
				}else{
					$best_d = $final_D;
				}
				print OUT "\t$best_d\t$best_j"; # print out
				if((split /\t/ ,$best_d)[0] ne "NA"){
					$D_align_n+=$abund_n;
					$vdj_n_stat{$seq_id}++ if(exists $vdj_n_stat{$seq_id} && $vdj_n_stat{$seq_id} == 2);
				}
				$pre_d = $_;
				%d_align_record = ();
				last;
			}
		}

		if((scalar keys %d_align_record) != 0) # the last records of D align file
		{
			my $best_d ;
			if($double_J)# re-alignment
			{
			$best_d = &best_align_d($seq , \%d_align_record , $v_end_pos , $j_start_pos , $strand_v); # select the best D align
			$best_d= join "\t" , (split /\s+/,$best_d)[1,2,3,4,6,7,8,9,10,11];
			}else{
				$best_d = $final_D;
			}
			print OUT "\t$best_d\t$best_j"; # print out
			if((split /\t/ ,$best_d)[0] ne "NA"){
				$D_align_n+=$abund_n;
				$vdj_n_stat{$seq_id}++ if(exists $vdj_n_stat{$seq_id} && $vdj_n_stat{$seq_id} == 2);
			}
			%d_align_record = ();
			$pre_d = "EOF";
		}
	}
	print OUT "\t$seq_id\t$seq\n";
}
close I;


#--------------------	alignment statistic
print "sequence_num: $fa_num\n";
print "V_alignment_rate: $V_align_n\t",$V_align_n/$fa_num*100,"\n";
print "D_alignment_rate: $D_align_n\t",$D_align_n/$fa_num*100,"\n";
print "J_alignment_rate: $J_align_n\t",$J_align_n/$fa_num*100,"\n";

my ($VJ_align_n ,$VDJ_align_n) = (0 , 0);
for(keys %vdj_n_stat){
	my $abund_n = (split /:/,$_)[1];
	$VDJ_align_n+=$abund_n if($vdj_n_stat{$_} == 3);
	$VJ_align_n+=$abund_n if($vdj_n_stat{$_} == 2 || $vdj_n_stat{$_} == 3);
}

print "VJ_alignment_rate: $VJ_align_n\t",$VJ_align_n/$fa_num*100,"\n";
print "VDJ_alignment_rate: $VDJ_align_n\t",$VDJ_align_n/$fa_num*100,"\n";


sub Modify_J_records
{
	my ($j_record, $str) = @_;
	return 1 if(scalar(keys %{$j_record}) == 1);
	my $max_score  = 0;
	for(keys %{$j_record})#get the best score
	{
		my @L = split;
		if($max_score < $L[-1]){
			$max_score = $L[-1];
		}
	}
	if($str eq "+")# 
	{
		my ($min_1,$min_2) = (1e4,1e4);
		for(keys %{$j_record})# get the minimal position
        	{	
			my @L = split;
			if($L[-1]> $max_score*0.8){
				($min_1,$min_2) = ($L[6],$L[7]) if($min_1 > $L[6]);
			}
		}
		for(keys %{$j_record})# delete incorrect J
        	{
                	my @L = split;
			if($L[6] > $min_2+50){
				delete($$j_record{$_});
			}
		}
	}
	else # minus chain
	{
		my ($max_1,$max_2) = (0,0);
	        for(keys %{$j_record})
        	{
                	my @L = split;
                	if($L[-1]> $max_score*0.8){
                        	($max_1,$max_2) = ($L[6],$L[7]) if($max_1 < $L[6]);
                	}
        	}
		for(keys %{$j_record})
		{
			my @L = split;
			if($L[7]+50 < $max_1){
                        	delete($$j_record{$_});
                	}
		}
	}
}



#-----------------------#
#	best_align_vj	#
#-----------------------#
sub best_align_vj
{
	my ($gene , $seq , $v_align , $edge_v , $strand_v) = @_;

	# V cannot find and $edge_v=NA
	if($gene eq "J" && $edge_v eq "NA")
	{
		if($strand_v eq "+"){$edge_v = 0;}
		else{$edge_v = length($seq)+1;}
	}

	my %best_align_all;# new align records
	my $max_score = -1000000;

	# 1. re-calculate the identity and the scorce, the new aligned start and end position

	for (keys %$v_align)
	{
		my $temp;
		next if((split /\s+/,$_)[7]-(split /\s+/,$_)[6] < $VJ_len_restrict);
		$temp = &re_calculate_score_identity($gene , $_ , $seq) if($gene eq "V"); # re-calculate subfunction
		if($gene eq "J") # re-calculate subfunction
		{
			$temp = &re_calculate_score_identity($gene , $_ , $seq , $edge_v , $strand_v);
		}
#	print "$temp\n" if($gene eq "J");

		next if((split /\s+/,$temp)[0] eq "NA");

		$max_score = (split /\s+/,$temp)[-1] if((split /\s+/,$temp)[-1] >= $max_score);
		push @{$best_align_all{(split /\s+/,$temp)[-1]}} , $temp;
	}
	

	# 2. select the best align record
	my $best_align;
	if(!exists $best_align_all{$max_score})# for J , the records were filter on the last step
	{
		return("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" , "NA" , "NA");
	}
	if(scalar @{$best_align_all{$max_score}} == 1){
		$best_align = $best_align_all{$max_score}->[0];
	}
	elsif(scalar @{$best_align_all{$max_score}} > 1)
	{
		$best_align = &select_best_align(@{$best_align_all{$max_score}});# have multiple max-score align records
	}
	# 3. find the V or J gene edge position and strand
	my ($edge_pos , $strand);
	if((split /\s+/,$best_align)[8]<(split /\s+/,$best_align)[9]){
		$edge_pos = (split /\s+/,$best_align)[7] if($gene eq "V");
		$edge_pos = (split /\s+/,$best_align)[6] if($gene eq "J");
		$strand = "+";
	}
	else{
		$edge_pos = (split /\s+/,$best_align)[6] if($gene eq "V");
		$edge_pos = (split /\s+/,$best_align)[7] if($gene eq "J");
		$strand = "-";
	}

	# 4. filter( identity && length) and return

	if($gene eq "V" && (split /\s+/,$best_align)[2] >= $v_identity_filter && (split /\s+/,$best_align)[3]>=$VJ_len_restrict)
	{
		return ($best_align , $edge_pos , $strand);
	}
	elsif($gene eq "J" && (split /\s+/,$best_align)[2] >= $j_identity_filter && (split /\s+/,$best_align)[3]>=$VJ_len_restrict)
	{
		return ($best_align , $edge_pos , $strand);
	}
	else
	{
		return ("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" , "NA" , "NA");
	}


}


#-----------------------#
#       best_align_d    #
#-----------------------#
sub best_align_d
{
	my ($seq , $d_align , $v_edge , $j_edge , $strand_v) = @_;
	if($v_edge eq "NA")# cannot find V align
	{
		return "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
	}
	if($j_edge eq "NA")# cannot find J align
	{
		$j_edge = 0 if($strand_v eq "-");
		$j_edge = length($seq)+1 if($strand_v eq "+");
	}

	my %best_align_all;# new align records
	my $max_score = 0;

	# 1. re-calculate the identity and the scorce, the new aligned start and end position
	for (keys %$d_align)
	{
		my $temp = &re_align_record_d($_ , $seq , $v_edge , $j_edge);# re-align
		
		if($temp ne "NA"){
			$max_score = (split /\s+/,$temp)[-1] if((split /\s+/,$temp)[-1] >= $max_score);
			push @{$best_align_all{(split /\s+/,$temp)[-1]}} , $temp;
		}
	}
	
	# 2. select the best align record
	my $best_align;
	if(!exists $best_align_all{$max_score})# for D , the records were filter on the last step
	{
		return("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
	}
	elsif(scalar @{$best_align_all{$max_score}} == 1){
		$best_align = $best_align_all{$max_score}->[0];
	}
	else{
		$best_align = &select_best_align_d(@{$best_align_all{$max_score}});# have multiple max-score align records
	}

	# 3. filter( identity ) and return
	if((split /\s+/,$best_align)[3] >= $D_len_restrict){
		return $best_align;
	}
	else{
		return ("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
	}
}


	


	#---------------------------------------#
	#	re-calculate score and identity	#
	#---------------------------------------#
sub re_calculate_score_identity
{
	
	my ($gene , $record , $seq , $edge_v , $strand_v) = @_ ;
	my $cdr3_p = (split /\./,(split /\s+/ , $record)[1])[-2];
	my @record_split = split /\s+/ , $record;


	my ($identity , $score);

	# 1. conduct sequences for calculating identity and score
	my @seq_base = split // , $seq;unshift(@seq_base , "N");
	my $ref_temp = $ref{$record_split[1]};

	# 1.1 for J, conduct the overlap region with V alignment
	my $J_ref_cut = 0;
	if($gene eq "J")
	{
		# store J align coordinate with ref
		
		my %J_align_check;
		my ($left_cor , $right_cor); # the coordinate of ref's edge

		if($record_split[8]<$record_split[9]){
			$left_cor = $record_split[6] - $record_split[8] + 1;
			$right_cor = $record_split[7] + length($ref_temp) - $record_split[9];
		}
		else{
			$left_cor = $record_split[6] - (length($ref_temp) - $record_split[8]);
			$right_cor = $record_split[7] + $record_split[9] - 1;
		}
		
		if(($strand_v eq "+" && $left_cor <= $edge_v) || ($strand_v eq "-" && $right_cor >= $edge_v))
		{

		if($record_split[8]<$record_split[9]){
			for(my $i=0 ; $i<=$record_split[7]-$record_split[6] ; $i++){$J_align_check{$record_split[6]+$i} = $record_split[8]+$i;}
		}else{
			for(my $i=0 ; $i<=$record_split[7]-$record_split[6] ; $i++){$J_align_check{$record_split[6]+$i} = $record_split[8]-$i;}
		}

		# conduct J align region
		if($strand_v eq "+")
		{
			if($record_split[6] <= $edge_v && $record_split[7] <= $edge_v) # all J in the V aligned region
			{
				return ("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" , "NA" , "NA");
			}
			elsif($record_split[6] <= $edge_v || $record_split[7] <= $edge_v)# part J in the V aligned region
			{
				my $l_c = $edge_v - $record_split[6]+1;
				my $n_r_c = $J_align_check{$edge_v};
				my @s_s = split // , (substr($seq , $record_split[6]-1 , $l_c));
				my @r_s;
				@r_s = split // , (substr($ref_temp , $record_split[8]-1 , $l_c)) if($record_split[8]<$record_split[9]);
				if($record_split[8]>$record_split[9])
				{
					my $t = substr($ref_temp , $n_r_c-1 , $l_c);
					$t =~ tr/ACGT/TGCA/;$t = reverse $t;
					@r_s = split // , $t;
				}

				my $mis_n = 0;
				for(my $i=0 ; $i<=$#r_s ; $i++){
					$mis_n++ if($s_s[$i] ne $r_s[$i]);
				}
				$record_split[3] -= $l_c;return ("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" , "NA" , "NA") if($record_split[3]<=6);
				$record_split[4] -= $mis_n;
				$record_split[6] = $edge_v+1;
				$record_split[8] = $J_align_check{$record_split[6]};
			}

			# change J ref sequence
			if($left_cor <= $edge_v)# change J ref sequence
			{
				if($record_split[8]<$record_split[9]){
					$ref_temp = substr($ref_temp , $edge_v-$left_cor+1);
					$J_ref_cut = $edge_v-$left_cor+1;
					$record_split[8] -= $J_ref_cut;
					$record_split[9] -= $J_ref_cut;
				}
				else{
					my $retain_len = length($ref_temp) - ($edge_v-$left_cor) - 1;
					return ("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA", "NA" , "NA") if($retain_len<=$cdr3_p);
					$ref_temp = substr($ref_temp , 0 , $retain_len);
				}
			}
			
		}
		else # for minus strand
		{
			if($record_split[6] >= $edge_v && $record_split[7] >= $edge_v) # all J in the V aligned region
			{
				return ("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA", "NA" , "NA");
			}
			elsif($record_split[6] >= $edge_v || $record_split[7] >= $edge_v)# part J in the V aligned region
			{
				my $l_c = $record_split[7] - $edge_v +1;
				my $n_r_c = $J_align_check{$edge_v};
				my @s_s = split // , (substr($seq , $edge_v-1 , $l_c));
				my @r_s;
				@r_s = split // , (substr($ref_temp , $n_r_c-1 , $l_c)) if($record_split[8]<$record_split[9]);
				if($record_split[8]>$record_split[9])
				{
					my $t = substr($ref_temp , $record_split[9]-1 , $l_c);
					$t =~ tr/ACGT/TGCA/; $t = reverse $t;
					@r_s = split // , $t;
				}
				my $mis_n = 0;
				for(my $i=0 ; $i<=$#r_s ; $i++){
					$mis_n++ if($s_s[$i] ne $r_s[$i]);
				}
				$record_split[3] -= $l_c;return ("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA", "NA" , "NA") if($record_split[3]<=6);
				$record_split[4] -= $mis_n;
				$record_split[7] = $edge_v-1;
				$record_split[9] = $J_align_check{$record_split[7]};
			}

			if($right_cor >= $edge_v)# change J ref sequence
			{
				if($record_split[8]>$record_split[9]){
					$ref_temp = substr($ref_temp , $right_cor-$edge_v+1);
					$J_ref_cut = $right_cor-$edge_v+1;
					$record_split[8] -= $J_ref_cut;
					$record_split[9] -= $J_ref_cut;
				}
				else{
					my $retain_len = length($ref_temp) - ($right_cor-$edge_v) - 1;
					return ("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA", "NA" , "NA") if($retain_len<=$cdr3_p);
					$ref_temp = substr($ref_temp , 0 , $retain_len);
				}
			}
		}
		}


	}
	
	my ($ref_s_p , $ref_e_p) = ($record_split[8] , $record_split[9]);
	if($record_split[8]>$record_split[9])# align to mimus strand
	{
		$ref_temp = reverse $ref_temp;$ref_temp=~tr/ACGT/TGCA/;# align to minus && ref sequence was translated to minus
		$ref_s_p = length($ref_temp)-$ref_s_p+1;# ref aligned position change
		$ref_e_p = length($ref_temp)-$ref_e_p+1;
	}
	my @ref_base = split // , $ref_temp;unshift @ref_base , "N";

	# 2. find the mismatch number on non-CDR3 region
	my $len_iden_score = $record_split[3];# the length for calculating identity and score
	my $mis_b_num = 0;# mismatch number
	if(($gene eq "V" && $record_split[8]<$record_split[9])||($gene eq "J" && $record_split[8]>$record_split[9]))
	{
#		$len_iden_score = $len_iden_score + $record_split[6] - 1;
		for (my $i=1; $i<$record_split[6] ; $i++) # 5' part
		{
			if($ref_s_p-$i >=1){
				$len_iden_score++;
				$mis_b_num++ if($ref_base[$ref_s_p-$i] ne $seq_base[$record_split[6]-$i]);
			}else{
				last; # the part exceeded V/J ref end will not be as a mismatch
			}
		}
		# change start alignment position 
		if($ref_s_p-$record_split[6]<0){
			if($gene eq "V"){
				$record_split[6] = $record_split[6]-$ref_s_p+1;
				$record_split[8] = 1;
			}
			else{
				$record_split[6] = $record_split[6]-$ref_s_p+1;
				$record_split[8] = length($ref_temp);
			}
		}
		else{
			if($gene eq "V"){
				$record_split[8] = $record_split[8]-$record_split[6]+1;
				$record_split[6] = 1;
			}else{
				$record_split[8] = $record_split[8]+$record_split[6]-1;
				$record_split[6] = 1;
			}
		}
	}
	elsif(($gene eq "V" && $record_split[8]>$record_split[9])||($gene eq "J" && $record_split[8]<$record_split[9]))
	{
		my $part_3_n = $#seq_base - $record_split[7];
#		$len_iden_score += $part_3_n;
		for (my $i=1 ; $i<=$part_3_n ; $i++)# 3'part
		{
			if($ref_e_p+$i <= $#ref_base){
				$len_iden_score++;
				$mis_b_num++ if($ref_base[$ref_e_p+$i] ne $seq_base[$record_split[7]+$i]);
			}else{
				last; #the part exceeded V/J ref end will not be as a mismatc
			}

#			$mis_b_num++ if($ref_e_p+$i > $#ref_base|| $ref_base[$ref_e_p+$i] ne $seq_base[$record_split[7]+$i]); # mismatch
		}
		# change start alignment position
		if($ref_e_p+$part_3_n > $#ref_base){
			if($gene eq "V"){
				$record_split[7] = $record_split[7]+$record_split[9]-1;
				$record_split[9] = 1;
			}else{
				$record_split[7] = $record_split[7]+(length($ref_temp)-$record_split[9]);
				$record_split[9] = length($ref_temp);
			}
		}
		else{
			if($gene eq "V"){
				$record_split[9] = $record_split[9]-($#seq_base-$record_split[7]);
				$record_split[7] = $#seq_base;
			}else{
				$record_split[9] = $record_split[9]+($#seq_base-$record_split[7]);
				$record_split[7] = $#seq_base;
			}
		}
	}
	$record_split[4] += $mis_b_num;
	$record_split[3] = $len_iden_score;


	# 3. re-alignment at CDR3 region
	my $record_new = join "\t" , @record_split;
	@record_split = &re_align_record($gene , $record_new , $v_mismatch , $seq , $J_ref_cut , $ref_temp) if($gene eq "V");
	@record_split = &re_align_record($gene , $record_new , $j_mismatch , $seq , $J_ref_cut , $ref_temp) if($gene eq "J");

	# 4. get the identity and score
	return "NA" if($record_split[3] <= 0);
	$identity = ($record_split[3]-$record_split[4])/$record_split[3]*100;$identity = sprintf("%0.2f",$identity);$record_split[2] = $identity;
	$score = $record_split[4]*$mismatch_s + ($record_split[3]-$record_split[4])*$match_s;$record_split[-1] = $score;
	
	$record_new = join "\t" , @record_split;

	return ($record_new);
}

	#-----------------------#
	#	re-align	#
	#-----------------------#

sub re_align_record
{
	my ($gene , $record , $mis_num , $seq , $J_ref_cut , $ref_temp) = @_ ;
	my $cdr3_p = (split /\./,(split /\s+/ , $record)[1])[-2];
	my @record_split = split /\s+/ , $record;


	my @seq_base = split // , $seq;unshift(@seq_base , "N");

	my ($ref_s_p , $ref_e_p) = ($record_split[8] , $record_split[9]);
	if($record_split[8]>$record_split[9])# align to mimus strand
	{
		$ref_s_p = length($ref_temp)-$ref_s_p+1;# ref aligned position change
		$ref_e_p = length($ref_temp)-$ref_e_p+1;
	}
	my @ref_base = split // , $ref_temp;unshift @ref_base , "N";


	# 4. re-assure the align edge
	my $cdr3_mis_raw = 0;
	my $cdr3_align_raw = 0;
	if($gene eq "V" && $record_split[8]<$record_split[9])# for V and forward strand
	{
		my $cdr3_mis_n = 0;
		if($ref_e_p >= $cdr3_p)# the align is more than CDR3 start 
		{
			my $flag = 0;
			$cdr3_align_raw = $ref_e_p-$cdr3_p+1;
			my $seq_cor_raw = $record_split[7];
			for (my $i=0 ; $i<=$#ref_base-$cdr3_p ; $i++)# 3'part
			{
				my $seq_index = $record_split[6]+$cdr3_p-$record_split[8]+$i;
				if( $seq_index >=1 && $seq_index<= $#seq_base &&  $ref_base[$cdr3_p+$i] ne $seq_base[$seq_index]){
					$cdr3_mis_n++;
					$cdr3_mis_raw++ if($seq_index<=$seq_cor_raw);
				}
				if($flag == 0 && ($seq_index <1 || $seq_index>$#seq_base ||$cdr3_mis_n == $mis_num+1))# 3'V-REGION reach $mis_num+1 mismatch
				{
					my $re_back = 0;
					if($i!=0 && $ref_base[$cdr3_p+$i-1] ne $seq_base[$seq_index-1]){
						my $sub_seq = substr($seq,$record_split[6]+$cdr3_p-$record_split[8]-1,$i);
						$sub_seq = reverse $sub_seq;$sub_seq=~tr/ACGT/TGCA/;
						my $sub_ref = substr($ref_temp , $cdr3_p-1 ,$i);
						$sub_ref=reverse$sub_ref;$sub_ref=~tr/ACGT/TGCA/;
						$re_back = &Remover_mis($sub_seq , $sub_ref);
					}

					$record_split[4] = $record_split[4]+$cdr3_mis_n-$re_back if($seq_index <1 || $seq_index>$#seq_base);
					$record_split[4] = $record_split[4]+$mis_num-$re_back if($cdr3_mis_n == $mis_num+1);
						
					$record_split[7] = $seq_index-1-$re_back;
					$record_split[9] = $cdr3_p+$i-1-$re_back;# new coordinate
					$record_split[3] = $record_split[3]+$i-$re_back;
					$flag = 1;
				}
			}
			if($flag == 0)# # 3'V-REGION don't reach $mis_num+1 mismatch
			{
				my $re_back = 0;
				if($ref_base[length($ref_temp)] ne $seq_base[$record_split[6]+$#ref_base-$record_split[8]]){
					my $sub_seq = substr($seq,$record_split[6]+$cdr3_p-$record_split[8]-1,$#ref_base-$cdr3_p+1);
					$sub_seq = reverse $sub_seq;$sub_seq=~tr/ACGT/TGCA/;
					my $sub_ref = substr($ref_temp , $cdr3_p-1 ,$#ref_base-$cdr3_p+1);
					$sub_ref=reverse$sub_ref;$sub_ref=~tr/ACGT/TGCA/;
					$re_back = &Remover_mis($sub_seq , $sub_ref);
				}
				$record_split[7] = $record_split[6] + $#ref_base-$record_split[8]-$re_back;
				$record_split[9] = $#ref_base-$re_back;
				$record_split[4] = $record_split[4]+$cdr3_mis_n-$re_back;
				$record_split[3] = $record_split[3]+$#ref_base-$cdr3_p+1-$re_back;
			}

		}
		else# the BLAST alignment exclude CDR3 region
		{
			my $distance = $cdr3_p-$record_split[9]-1;
			my $mismatch_n = 0;
			for(my $i=1 ; $i<=$distance ; $i++){
				$mismatch_n++ if(($i+$record_split[7])>$#seq_base || ($record_split[9]+$i)>$#ref_base || $seq_base[$i+$record_split[7]] ne $ref_base[$record_split[9]+$i]);	
			}
			$record_split[3] = $record_split[3]+$distance;
			$record_split[4] = $record_split[4]+$mismatch_n;
		}

	}

	elsif($gene eq "V" && $record_split[8]>$record_split[9]) # for V and minus strand
	{

		my $cdr3_mis_n = 0;
		if($record_split[8] >= $cdr3_p)# the align is more than CDR3 start 
		{
			my $flag = 0;
			$cdr3_align_raw = $record_split[8]-$cdr3_p+1;
			my $seq_cor_raw = $record_split[6];
			for (my $i=0 ; $i<=$#ref_base-$cdr3_p ; $i++)# 3'part
			{
				my $seq_index = $record_split[7]-$cdr3_p+$record_split[9]-$i;
				if($seq_index >=1 && $seq_index<= $#seq_base && $ref_base[$#ref_base-$cdr3_p-$i+1] ne $seq_base[$seq_index]){
					$cdr3_mis_n++;
					$cdr3_mis_raw++ if($seq_index>=$seq_cor_raw);
				}
				if($flag == 0 && ($seq_index <1 || $seq_index>$#seq_base || $cdr3_mis_n == $mis_num+1))# 3'V-REGION reach $mis_num+1 mismatch
				{
					my $re_back = 0;
					if($i!=0 && $ref_base[$#ref_base-$cdr3_p-$i+2] ne $seq_base[$seq_index+1]){
						my $sub_seq = substr($seq,$seq_index,$i);
						my $sub_ref = substr($ref_temp,$#ref_base-$cdr3_p-$i+1,$i);
						$re_back = &Remover_mis($sub_seq , $sub_ref);
					}
					$record_split[4] = $record_split[4]+$cdr3_mis_n-$re_back if($seq_index <1 || $seq_index>$#seq_base);
					$record_split[4] = $record_split[4]+$mis_num-$re_back if($cdr3_mis_n == $mis_num+1);

					$record_split[6] = $seq_index+1+$re_back;# new coordinate
					$record_split[8] = $cdr3_p+$i-1-$re_back;# new coordinate
					$record_split[3] = $record_split[3]+$i-$re_back;
					$flag = 1;
				}
			}
			if($flag == 0)# # 3'V-REGION don't reach $mis_num+1 mismatch
			{
				my $re_back = 0;
				if($ref_base[1] ne $seq_base[$record_split[7]+$record_split[9]-$#ref_base]){
					my $sub_seq = substr($seq,$record_split[7]+$record_split[9]-$#ref_base-1,$#ref_base-$cdr3_p+1);
					my $sub_ref = substr($ref_temp,0,$#ref_base-$cdr3_p+1);
					$re_back = &Remover_mis($sub_seq , $sub_ref);
				}
				$record_split[6] = $record_split[7] - ($#ref_base - $record_split[9])+$re_back;
				$record_split[8] = $#ref_base-$re_back;
				$record_split[4] = $record_split[4]+$cdr3_mis_n-$re_back;
				$record_split[3] = $record_split[3]+$#ref_base-$cdr3_p+1-$re_back;
			}
		}
		else# the BLAST alignment exclude CDR3 region
		{
			my $distance = $cdr3_p - $record_split[8]-1;
			my $mismatch_n = 0;
			for(my $i=1 ; $i<=$distance ; $i++){
				$mismatch_n++ if(($record_split[6]-$i)<=0 || $ref_s_p-$i<=0 || $seq_base[$record_split[6]-$i] ne $ref_base[$ref_s_p-$i]);
			}
			$record_split[3] = $record_split[3] + $distance;
			$record_split[4] = $record_split[4] + $mismatch_n;
		}
	}
	elsif($gene eq "J" && $record_split[8]<$record_split[9])# for J and forward strand
	{
	
		my $cdr3_mis_n = 0;
		$cdr3_p -= $J_ref_cut;
		if($ref_s_p <= $cdr3_p)# the align is more than CDR3 start 
		{
			my $flag = 0;
			$cdr3_align_raw = $cdr3_p-$ref_s_p+1;
			my $seq_cor_raw = $record_split[6];
			for (my $i=0 ; $i<$cdr3_p ; $i++)# 5'part
			{
				my $seq_index = $record_split[7]-($record_split[9]-$cdr3_p)-$i;
				if($seq_index >=1 && $seq_index<= $#seq_base && $ref_base[$cdr3_p-$i] ne $seq_base[$seq_index]){
					$cdr3_mis_n++;
					$cdr3_mis_raw++ if($seq_index>=$seq_cor_raw);
				}
				if($flag == 0 && ($seq_index <1 || $seq_index>$#seq_base || $cdr3_mis_n == $mis_num+1))# 5'J-REGION reach $mis_num+1 mismatch
				{
					my $re_back = 0;
					if($i!=0 && $ref_base[$cdr3_p-$i+1] ne $seq_base[$seq_index+1]){
						my $sub_seq = substr($seq,$seq_index,$i);
						my $sub_ref = substr($ref_temp,$cdr3_p-$i,$i);
						$re_back = &Remover_mis($sub_seq , $sub_ref);
					}
					$record_split[4] = $record_split[4]+$cdr3_mis_n-$re_back if($seq_index <1 || $seq_index>$#seq_base);
					$record_split[4] = $record_split[4]+$mis_num-$re_back if($cdr3_mis_n == $mis_num+1);

					$record_split[6] = $seq_index+1+$re_back;# new coordinate
					$record_split[8] = $cdr3_p-$i+1+$re_back;# new coordinate
					$record_split[3] = $record_split[3]+$i-$re_back;
					$flag = 1;
				}
			}
			if($flag == 0)# # 5'J-REGION don't reach $mis_num+1 mismatch
			{
				my $re_back = 0;
				if($ref_base[1] ne $seq_base[$record_split[7]-$record_split[9]+1]){
					my $sub_seq = substr($seq,$record_split[7]-$record_split[9],$cdr3_p);
					my $sub_ref = substr($ref_temp,0,$cdr3_p);
					$re_back = &Remover_mis($sub_seq , $sub_ref);
				}
				$record_split[6] = $record_split[7]-$record_split[9]+1+$re_back;
				$record_split[8] = 1+$re_back;
				$record_split[4] = $record_split[4]+$cdr3_mis_n-$re_back;
				$record_split[3] = $record_split[3]+$cdr3_p-$re_back;
			}
		}
		else# the BLAST alignment exclude CDR3 region
		{
			my $distance = $record_split[8]-$cdr3_p-1;
			my $mismatch_n = 0;
			for (my $i=0 ; $i<=$distance ; $i++){
				$mismatch_n++ if(($record_split[6]-$i)<=0 || $ref_s_p-$i<=0 || $seq_base[$record_split[6]-$i] ne $ref_base[$ref_s_p-$i]);
			}
			$record_split[3] = $record_split[3]+$distance;
			$record_split[4] = $record_split[4]+$mismatch_n;
		}
	}
	elsif($gene eq "J" && $record_split[8]>$record_split[9])# for J and minus strand
	{

		my $cdr3_mis_n = 0;
		$cdr3_p -= $J_ref_cut;
		if($record_split[9] <= $cdr3_p)# the align is more than CDR3 start 
		{
			my $flag = 0;
			$cdr3_align_raw = $cdr3_p-$record_split[9]+1;
			my $seq_cor_raw = $record_split[7];
			for (my $i=0 ; $i<$cdr3_p ; $i++)# 5'part
			{
				my $seq_index = $record_split[6]+$record_split[8]-$cdr3_p+$i;
				if($seq_index >=1 && $seq_index<= $#seq_base && $ref_base[$#ref_base-$cdr3_p+$i+1] ne $seq_base[$seq_index]){
					$cdr3_mis_n++;
					$cdr3_mis_raw++ if($seq_index<=$seq_cor_raw);
				}
				if($flag == 0 && ($seq_index <1 || $seq_index>$#seq_base || $cdr3_mis_n == $mis_num+1))# 5'J-REGION reach $mis_num+1 mismatch
				{
					my $re_back = 0;
					if($i!=0 && $ref_base[$#ref_base-$cdr3_p+$i] ne $seq_base[$seq_index-1]){
						my $sub_seq = substr($seq,$record_split[6]+$record_split[8]-$cdr3_p-1,$i);
						$sub_seq=reverse $sub_seq;$sub_seq=~tr/ACGT/TGCA/;
						my $sub_ref = substr($ref_temp,$#ref_base-$cdr3_p,$i);
						$sub_ref=reverse $sub_ref;$sub_ref=~tr/ACGT/TGCA/;
						$re_back = &Remover_mis($sub_seq , $sub_ref);
					}
					$record_split[4] = $record_split[4]+$cdr3_mis_n-$re_back if($seq_index <1 || $seq_index>$#seq_base);
					$record_split[4] = $record_split[4]+$mis_num-$re_back if($cdr3_mis_n == $mis_num+1);

					$record_split[7] = $seq_index-1-$re_back;# new coordinate
					$record_split[9] = $cdr3_p-$i+1+$re_back;# new coordinate
					$record_split[3] = $record_split[3]+$i-$re_back;
					$flag = 1;
				}
			}
			if($flag == 0)# # 5'J-REGION don't reach $mis_num+1 mismatch
			{
				my $re_back = 0;
				if($ref_base[length($ref_temp)] ne $seq_base[$record_split[6]+$record_split[8]-1]){
					my $sub_seq = substr($seq,$record_split[6]+$record_split[8]-$cdr3_p-1,$cdr3_p);
					$sub_seq=reverse $sub_seq;$sub_seq=~tr/ACGT/TGCA/;
					my $sub_ref = substr($ref_temp,$#ref_base-$cdr3_p,$cdr3_p);
					$sub_ref=reverse $sub_ref;$sub_ref=~tr/ACGT/TGCA/;
					$re_back = &Remover_mis($sub_seq , $sub_ref);
				}
				$record_split[7] = $record_split[6]+$record_split[8]-1-$re_back;
				$record_split[9] = 1+$re_back;
				$record_split[4] = $record_split[4]+$cdr3_mis_n-$re_back;
				$record_split[3] = $record_split[3]+$cdr3_p-$re_back;
			}
		}
		else# the BLAST alignment exclude CDR3 region
		{
			my $distance = $record_split[9] - $cdr3_p - 1;
			my $mismatch_n = 0;
			for(my $i=1 ; $i<=$distance ; $i++){
				$mismatch_n++ if(($record_split[7]+$i)>$#seq_base || $ref_e_p+$i>$#ref_base || $seq_base[$record_split[7]+$i] ne $ref_base[$ref_e_p+$i]);
			}
			$record_split[3] = $record_split[3] + $distance;
			$record_split[4] = $record_split[4] + $mismatch_n;
		}
	}
	$record_split[4] = $record_split[4]-$cdr3_mis_raw;
	$record_split[3] = $record_split[3]-$cdr3_align_raw;
	if($gene eq "J" && $J_ref_cut)
	{
		$record_split[8] += $J_ref_cut;
		$record_split[9] += $J_ref_cut;
	}
	return (@record_split);
}



#------------------------------#
#	align edge conduct     #
#------------------------------#
sub Remover_mis
{
	my ($s1,$s2) = @_;
	my $re_back = 0;
	my @s_1 = split //,$s1;
	my @s_2 = split //,$s2;
	for(my $i=0;$i<=$#s_1;$i++)
	{
		last if($s_1[$i] eq $s_2[$i]);
		$re_back++ if($s_1[$i] ne $s_2[$i]);
	}
	return $re_back;
}

#-------------------------------#
#	re_align_record_d	#
#-------------------------------#
sub re_align_record_d
{
	my ($record , $seq , $edge_v , $edge_j) = @_;
	my @record_split = split /\s+/ , $record;
	($edge_v , $edge_j) = sort {$a<=>$b} ($edge_v , $edge_j);

	my ($identity , $score);
	# 1. conduct sequences for calculating identity and score
	my @seq_base = split // , $seq;unshift(@seq_base , "N");
	my $ref_temp = $ref{$record_split[1]};
	my ($ref_s_p , $ref_e_p) = ($record_split[8] , $record_split[9]);
	if($record_split[8]>$record_split[9])# align to mimus strand
	{
		$ref_temp = reverse $ref_temp;$ref_temp=~tr/ACGT/TGCA/;# align to minus && ref sequence was translated to minus
		$ref_s_p = length($ref_temp)-$ref_s_p+1;# ref aligned position change
		$ref_e_p = length($ref_temp)-$ref_e_p+1;
	}
	my @ref_base = split // , $ref_temp;unshift @ref_base , "N";
	
	# 2. re-assure the align start and end position
	
	if($record_split[6]+$#ref_base-$ref_e_p<= $edge_v || $record_split[6]-$ref_s_p+1 >= $edge_j) # the all D seq align to V region or J region
	{
		return "NA";
	}

	my ($new_s_q , $new_e_q, $new_s_r , $new_e_r);
	if($record_split[6]-$ref_s_p+1 <= $edge_v) # the start ref position is more than $edge_v
	{
		$new_s_q = $edge_v+1;
		$new_s_r = $ref_s_p - ($record_split[6]-$edge_v)+1;
	}
	else
	{
		$new_s_q = $record_split[6]-$ref_s_p+1;
		$new_s_r = 1;
	}

	if($record_split[7]+$#ref_base-$ref_e_p >= $edge_j) # the end ref position is more than $edge_j
	{
		$new_e_q = $edge_j - 1;
		$new_e_r = $ref_e_p + $edge_j-$record_split[7]-1;
	}
	else
	{
		$new_e_q = $record_split[7]+$#ref_base-$ref_e_p;
		$new_e_r = $#ref_base;
	}

	# 3. re-align and select the best align result according mismatch number
	my %seq_ref_check;
	my $k = 0;
	for ($new_s_q..$new_e_q) # the base_seq 's cordinate -> the ref 's cordinate
	{
		$seq_ref_check{$_} = $new_s_r+$k;
		$k++;
	}

	my %re_align_all;
	for (my $i=$new_s_q ; $i<=$new_e_q ; $i++)
	{
		if($i==$new_s_q || $seq_base[$i] ne $ref_base[$seq_ref_check{$i}])# the start position for calculating score
		{
			# the start position
			my ($seq_start,$j);
			if($i==$new_s_q){$j=$i;}
			else{$j=$i+1;}
			$seq_start = $j;
			
			my $mis_num = 0;
			for (;$j<=$new_e_q; $j++)
			{
				$mis_num++ if($seq_base[$j] ne $ref_base[$seq_ref_check{$j}]);
				if($mis_num == $d_mismatch+1) # reach the D mismatch restrict
				{
					my $len = $j-$seq_start;
					my $score = $d_mismatch*$mismatch_s+($len-$d_mismatch)*$match_s;
					$re_align_all{"$seq_start:".($j-1).":$len:$d_mismatch"} = $score if($score != 0);
					last;
				}
			}
			if($mis_num <= $d_mismatch)
			{
				my $len = $j-$seq_start;
				my $score = $mis_num*$mismatch_s+($len-$mis_num)*$match_s;
				$re_align_all{"$seq_start:".($j-1).":$len:$mis_num"} = $score if($score != 0);
			}
		}
	}

	# 4. select the best align, filter and return
	return "NA" if((scalar keys %re_align_all) == 0);

	my $best = (sort {$re_align_all{$b}<=>$re_align_all{$a}} keys %re_align_all)[0];
	($record_split[6] , $record_split[7] , $record_split[3] , $record_split[4]) = split /:/,$best;

	$record_split[2] = (1 - $record_split[4]/$record_split[3])*100;$record_split[2] = sprintf("%0.2f",$record_split[2]);

	return "NA" if($record_split[3] < $D_len_restrict);# use align len to filter
	
	my ($ref_start , $ref_end) = ($seq_ref_check{$record_split[6]} , $seq_ref_check{$record_split[7]});
	if($record_split[8] > $record_split[9])# minus align and change cordinate
	{
		$record_split[8] = $#ref_base-$ref_start+1;
		$record_split[9] = $#ref_base-$ref_end+1;
	}
	else
	{
		($record_split[8] , $record_split[9]) = ($ref_start , $ref_end);
	}
	$record_split[-1] = $re_align_all{$best};
	$best = join "\t" , @record_split;
}


	#-----------------------#
	#	best_align_all	#
	#-----------------------#
	# for same score records, first select the minimum unalignment length of CDR3 edge , then rand the leftover records

sub select_best_align
{
	my $best_align;
	my @b_a_a = @_;
	my $max_id = 1000;
	my %b_a_a_d;
	for (@b_a_a)
	{
		my $ref_max = (split /\s+/,$_)[8]>(split /\s+/,$_)[9]?(split /\s+/,$_)[8]:(split /\s+/,$_)[9];
		my $unmap_len = length($ref{(split /\s+/,$_)[1]})-$ref_max;

		push @{$b_a_a_d{$unmap_len}} , $_;
		$max_id = $unmap_len if($unmap_len <= $max_id);
	}
	$best_align = $b_a_a_d{$max_id}->[int rand(scalar @{$b_a_a_d{$max_id}})];
	return $best_align;
}
        #-----------------------#
        #       best_align_all  #
        #-----------------------#
        # for same score records, first select the max identity records, then rand the leftover records

sub select_best_align_d
{
        my $best_align;
        my @b_a_a = @_;
        my $max_id = 0;
        my %b_a_a_d;
        for (@b_a_a)
        {
                push @{$b_a_a_d{(split /\s+/,$_)[2]}} , $_;
                $max_id = (split /\s+/,$_)[2] if((split /\s+/,$_)[2] >= $max_id);
        }
        $best_align = $b_a_a_d{$max_id}->[int rand(scalar @{$b_a_a_d{$max_id}})];
        return $best_align;
}

