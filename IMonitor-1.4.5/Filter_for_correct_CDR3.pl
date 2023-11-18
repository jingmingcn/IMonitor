#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);

die "perl $0 <ref>  <in.gz> <fq.gz>  <id.backup> <out.gz> \n" unless(@ARGV==5);

my($ref , $in , $fq  , $id_backup, $out) = @ARGV;

#my $p_mis_n=2;
#my $primer_l_restrict;
#unless(defined $primer_l_restrict){$primer_l_restrict=10;}

open O , "|gzip >$out" or die;

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


my %id_change;
open I, "gzip -dc $id_backup|" or die;
while(<I>)
{
	chomp;
	my @line = split;
	$id_change{$line[0]} = $line[1];
}
close I;

if($fq=~/\.gz$/){
	open Q, "gzip -dc $fq|" or die;
}else{
	open Q, "$fq" or die;
}
#%id_change=();

$in!~/\.gz$/ and (open I, $in or die "Failed in opening $in !!\n");
$in=~/\.gz$/ and (open I, "gzip -dc $in|" or die "Failed in opening $in !!\n");

#open D, ">$discard_f" or die;


#print D "##Reason_of_Discarding: 1\tV_or_J_gene_missed\n##Reason_of_Discarding: 3\tAligned_to_conflicted_strands\n##Reason_of_Discarding: 5\tFind primer length less than 10+bps or exists more than 2 mismatch\n";

my ($vj_f , $vj_conf, $cdr3_err) = (0,0,0);

while (<I>){
	chomp;
	my $raw = $_;
	my @F=split /\s+/;
	
	# 1. filter the sequences
	if(($F[0] eq "NA")||($F[20] eq "NA"))# filter the records without V or J gene
	{
#		print D "1\t$raw\n";
		$vj_f++;
		next;
	}
	elsif(($F[6]-$F[7])*($F[26]-$F[27])<0)# filter the records with V and J alignment conflict
	{
#		print D "2\t$raw\n";
		$vj_conf++;
		next;
	}

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

	if(($F[6]-$F[7])>0) #change the minus to forward strand
	{
		@F=&trans_complementary_str(@F);
		$qual = reverse $qual;
	}


	# 3. find the CDR3 
	my $cdr3_ref_start = (split /\./,$F[0])[4];
	my $cdr3_ref_end = (split /\./,$F[20])[4];
	my $cdr3_start = $F[5]-($F[7]-$cdr3_ref_start);
	my $cdr3_end = $F[24] + ($cdr3_ref_end-$F[26]);
	
	my $cdr3length=$cdr3_end-$cdr3_start+1;
	if($cdr3length<=0)# CDR3 length error
	{	
		$cdr3_err++;
#		print D "4\t$raw\n";
		next;
	}
	
	my $cdr3 = substr($F[31],$cdr3_start-1,$cdr3length);
	my $cdr3_qual = substr($qual,$cdr3_start-1,$cdr3length);


	print O join "\t", @F;
	print O "\t$cdr3\t$cdr3_qual\t0\n";

	# 
}
close I;


print "without_V_or_J: $vj_f\n";
print "V_and_J_strand_conflict: $vj_conf\n";
print "CDR3_length_error: $cdr3_err\n";


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

