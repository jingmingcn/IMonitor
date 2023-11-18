#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);

die "perl $0 <ref> <in.gz> <fq.gz> <id.backup>\n" unless(@ARGV==4);

my($ref , $in , $fq , $id_backup) = @ARGV;

my $p_mis_n=2;
my $primer_l_restrict;
unless(defined $primer_l_restrict){$primer_l_restrict=10;}

#open O , "|gzip >$out" or die;

my %ref;
my $name;
#Get reference information: $ref{reference name}= sequence of reference
open R, $ref or die "Failed in opening $ref !!\n";
while (<R>){
        chomp;
        /^>([^\n]+)/ and $name=$1 and $ref{$name}="";
        !/^>/ and $ref{$name}.=$_;
}
close R;

#print Dumper(\%primer_tem_check);

my %id_change;
open I, "gzip -dc $id_backup|" or die;
while(<I>)
{
	chomp;
	my @line = split;
	$id_change{$line[0]} = $line[1];
}
close I;

open Q, "$fq" or die;
#%id_change=();

$in!~/\.gz$/ and (open I, $in or die "Failed in opening $in !!\n");
$in=~/\.gz$/ and (open I, "gzip -dc $in|" or die "Failed in opening $in !!\n");

#open D, ">$discard_f" or die;


#print D "##Reason_of_Discarding: 1\tV_or_J_gene_missed\n##Reason_of_Discarding: 3\tAligned_to_conflicted_strands\n##Reason_of_Discarding: 5\tFind primer length less than 10+bps or exists more than 2 mismatch\n";

my ($vj_f , $vj_conf) = (0,0);

while (<I>){
	chomp;
	my $raw = $_;
	my @F_raw = split /\s+/,$_;
	my @F= @F_raw[19..48,0,17];
	my $strand_f = $F_raw[16];
	
	# 1. filter the sequences
#	if(($F[0] eq "NA")||($F[20] eq "NA"))# filter the records without V or J gene
#	{
#		print D "1\t$raw\n";
#		$vj_f++;
#		next;
#	}
#	elsif(($F[6]-$F[7])*($F[26]-$F[27])<0)# filter the records with V and J alignment conflict
#	{
#		print D "2\t$raw\n";
#		$vj_conf++;
#		next;
#	}

	# 2. conduct the sequences
	my $qual;
	while(my $id = <Q>)
	{
		chomp($id);
		$id = (split /\s+/,$id)[0];
		<Q>;<Q>;
		chomp($qual=<Q>);
		last if(exists $id_change{$id} and $id_change{$id} eq $F[30]);
	}

	if($strand_f eq "Minus") #change the minus to forward strand
	{
#		@F=&trans_complementary_str(@F);
		$qual = reverse $qual;
	}

	my $v_add_l = 0;


	# for 5'RACE, add the V not sequenced and change the minus to forward strand
	
		my $new_seq = $F[31];
#		if($F[6]-$F[4]>0) # add the part V sequence
#		{
#			my $missed_V_primer_seq=substr($ref{$F[0]},0,$F[6]-$F[4]);
#			$new_seq = $missed_V_primer_seq.$new_seq;
#			$v_add_l = length($missed_V_primer_seq);
#		}
		$qual = ("h" x $v_add_l).$qual;
		print join "\t", @F;
		print "\t$new_seq\t$qual\t$v_add_l\n";
	

	# 
}
close I;


#print "without_V_or_J: $vj_f\n";
#print "V_and_J_strand_conflict: $vj_conf\n";


sub trans_complementary_str{
	my @F = @_;
	my $tmp;
	my $q_len=length($F[31]);

	$F[31]=reverse $F[31];
	$F[31]=~tr/ATCGatcg/TAGCtagc/;
	$F[4]=$q_len-$F[4]+1;
	$F[5]=$q_len-$F[5]+1;
	$F[24]=$q_len-$F[24]+1;
	$F[25]=$q_len-$F[25]+1;
	if($F[10] ne "NA"){
		$F[14]=$q_len-$F[14]+1;
		$F[15]=$q_len-$F[15]+1;
		$tmp=$F[14];$F[14]=$F[15];$F[15]=$tmp;
		$tmp=$F[16];$F[16]=$F[17];$F[17]=$tmp;
	}
	
	$tmp=$F[4];$F[4]=$F[5];$F[5]=$tmp;
	$tmp=$F[24];$F[24]=$F[25];$F[25]=$tmp;
	$tmp=$F[6];$F[6]=$F[7];$F[7]=$tmp;
	$tmp=$F[26];$F[26]=$F[27];$F[27]=$tmp;

	return(@F);
}

#-----------------------#
#	cut primer seq	#
#-----------------------#
#
# 1. don't cut primer: 1). V/J don't have align result(NA); 
#                      2). V and J align to different strands;
#		       3). error align: forward strand (J at the left of V)/ minus strand(V at the left of J)
# 2. return "NA	NA ...": V or J aligned region all during the primer region

	#-------------------------------------------
	#	2 sequence alignment with mismatch
	#------------------------------
sub primer_align
{
	my ($s1 , $s2 , $mis_n) = @_;
	return 1 if($s1 eq $s2);
	my $mismatch = 0;
	my @s_1 = split //,$s1;
	my @s_2 = split //,$s2;
	for (my $i=0 ; $i<=$#s_1 ; $i++)
	{
		$mismatch++ if($s_1[$i] ne $s_2[$i]);
		return 0 if($mismatch >$mis_n);
	}
	return 1;
}
