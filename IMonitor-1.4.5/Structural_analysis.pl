#!/usr/bin/perl -w
use strict;
use Data::Dumper;

=head1 Usage

	Structural_analysis.pl  -i <> -o <> -n <> -r <> [parameters] 

=head1 Parameters

        -i	input *.SingleLine file
        -o	output directory path [current directory]
        -n	sample name [basename of -b arguemtn]
	-r	reference fa file
	-d	logical value, if add -d, need to analyze D gene
	-f1	CDR3(nucleotide) abundance filter [0]
	-f2	Clonotype(nucleotide, ful-length) abundance filter [0]
	
	If it need to find the primer, then use those parameters as follow: -u [-p -m]
	-u	Primer sequence and template sequence check file
	-p      the number of primers length restrict, sequences less than the length will be discarded [10]
	-m      mismatch of found primer sequence [2]

	-h	print help information

=head1 Function
	1. filter the records without V or J gene
	2. filter the records with V and J align strand conflict
	3. filter the CDR3 length <= 0
	4. optionally, find the multiplex primers and filter the sequence witout primer
	5. change the minus strand of records to plus strand
	6. Find CDR3 region
	7. Translate nucleotide sequence into animo acid 
	8. add a flag: "out-of-frame","in-frame","non-function".(only one of V or J is a peseudogene or ORF, then give a "non-function")
	9. optionally, sequence with low abundance can be filtered: -f1, -f2


=cut

use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
my ($in,$out_dir,$ref,$ind,$inb,$sample_name,$primer_l_restrict,$p_mis_n ,$Help);
my ($filter_cdr3,$filter_clone) = (0,0);

GetOptions(
        "i:s" => \$in,
        "o:s" => \$out_dir,
        "n:s" => \$sample_name,
        "p:s" => \$primer_l_restrict,
	"m:i" => \$p_mis_n,
        "r:s" => \$ref,
	"d" => \$ind,
	"u:s" => \$inb,
	"f1:i" => \$filter_cdr3,
	"f2:i" => \$filter_clone,
        "help!" => \$Help,
);

$out_dir||=".";

`mkdir -p $out_dir` unless(-d $out_dir);
$p_mis_n=2 unless(defined($p_mis_n));

unless(defined $primer_l_restrict){$primer_l_restrict=10;}

die `pod2text $0` if ($Help or !$in or !$out_dir or !$sample_name or !$ref);


my %codon2aa; # amino acid codon
        $codon2aa{"TTT"}= "F"; $codon2aa{"TTC"}= "F"; $codon2aa{"TTA"}="L";
        $codon2aa{"TTG"}= "L"; $codon2aa{"TCT"}= "S"; $codon2aa{"TCC"}= "S";
        $codon2aa{"TCA"}= "S"; $codon2aa{"TCG"}= "S"; $codon2aa{"TAT"}= "Y";
        $codon2aa{"TAC"}= "Y"; $codon2aa{"TGT"}= "C"; $codon2aa{"TGC"}= "C";
        $codon2aa{"CTT"}= "L"; $codon2aa{"CTC"}= "L"; $codon2aa{"CTA"}= "L";
        $codon2aa{"CTG"}= "L"; $codon2aa{"CCT"}= "P"; $codon2aa{"CCC"}= "P";
        $codon2aa{"CCA"}= "P"; $codon2aa{"CCG"}= "P"; $codon2aa{"CAT"}= "H";
        $codon2aa{"CAC"}= "H"; $codon2aa{"CAA"}= "Q"; $codon2aa{"CAG"}= "Q";
        $codon2aa{"CGT"}= "R"; $codon2aa{"CGC"}= "R"; $codon2aa{"CGA"}= "R";
        $codon2aa{"CGG"}= "R"; $codon2aa{"ATT"}= "I"; $codon2aa{"ATC"}= "I";
        $codon2aa{"ATA"}= "I"; $codon2aa{"ATG"}= "M"; $codon2aa{"ACT"}= "T";
        $codon2aa{"ACC"}= "T"; $codon2aa{"ACA"}= "T"; $codon2aa{"ACG"}= "T";
        $codon2aa{"AAT"}= "N"; $codon2aa{"AAC"}= "N"; $codon2aa{"AAA"}= "K";
        $codon2aa{"AAG"}= "K"; $codon2aa{"AGT"}= "S"; $codon2aa{"AGC"}= "S";
        $codon2aa{"AGA"}= "R"; $codon2aa{"AGG"}= "R"; $codon2aa{"GTT"}= "V";
        $codon2aa{"GTC"}= "V"; $codon2aa{"GTA"}= "V"; $codon2aa{"GTG"}= "V";
        $codon2aa{"GCT"}= "A"; $codon2aa{"GCC"}= "A"; $codon2aa{"GCA"}= "A";
        $codon2aa{"GCG"}= "A"; $codon2aa{"GAT"}= "D"; $codon2aa{"GAC"}= "D";
        $codon2aa{"GAA"}= "E"; $codon2aa{"GAG"}= "E"; $codon2aa{"GGT"}= "G";
        $codon2aa{"GGC"}= "G"; $codon2aa{"GGA"}= "G"; $codon2aa{"GGG"}= "G";
        $codon2aa{"TGG"}= "W";
        $codon2aa{"TAG"}= "*";
        $codon2aa{"TGA"}= "*";
        $codon2aa{"TAA"}= "*";


#output file of detailed information
my $out="$out_dir/$sample_name.structure";
#output file of record without reasonable alignment(missing one or two of VDJ genes, or conflicted alignment strand)
my $err="$out_dir/$sample_name.discarded";
#

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

my %primer_tem_check;
if(defined($inb))
{
	open I, "$inb" or die;
	while(<I>)
	{
		chomp;
		my @line = split;
		$primer_tem_check{$line[0]} = $line[3];
	}
	close I;
}


#my ($total_seq,$vdj_seq,$orf_seq,$D_inversion)=(0,0,0,0);
my %discard;
my %discard_header=(
	1=>"V_or_J_gene_missed",
	2=>"Aligned_to_conflicted_strands",
	3=>"",
	4=>"CDR3_length_<=0bp",
);


$in!~/\.gz$/ and (open I, $in or die "Failed in opening $in !!\n");
$in=~/\.gz$/ and (open I, "gunzip -c $in|" or die "Failed in opening $in !!\n");
open O, "|gzip -c >$out.gz" or die "Failed in creating $out.gz !!\n";
open E, "|gzip -c >$err.gz" or die "Failed in creating $err.gz !!\n";
print O "#ID\tfunction\tV_ref\tD_ref\tJ_ref\tCDR3_start\tCDR3_end\tCDR3(dna)\tCDR3(aa)\t3'V_del\t5'D_del\t3'D_del\t5'J_del\tVD_ins\tDJ_ins\tVJ_ins\tstrand\tsequence\tamino_acid\talignment_record\n";
open T, ">$out_dir/$sample_name.structure.tmp.o" or die;


for my $err_ind(1..4){
	print E "##Reason_of_Discarding: $err_ind\t$discard_header{$err_ind}\n";
}

#------------------------	read *.Single.gz file	-----------------------
my ($total_seq, $primer_filter , $without_vj , $vj_conflict , $cdr3_err) = (0,0,0,0,0);
my %cdr3_abund;
my %clone_abund;

while (<I>){
        chomp;
        my $ori_line=$_;
#$.==10000 and print "$vdj_seq\n"and exit;
#Fields: 1:Subject id, 2:% identity, 3:alignment length, 4:mismatches, 5:q. start, 6:q. end, 7:s. start, 8:s. end, 9:e-value, 10:bit score ... 31:query id, 32:query sequence
#1:IGHV3-23*03/IGHV3-23*05        2:98.04        3:51        4:1        5:65        6:115        7:51        8:1        9:9e-23        10:93.7        11:IGHD6-13*01        12:93.75        13:16        14:1        15:3816:53        17:21        18:6        19:0.005        20:24.3        21:IGHJ4*03        22:94.12        23:34        24:2        25:1        26:34        27:41        28:8        29:3e-11        30:52.0        31:cp1_85_115        32:GACGGTGACCAGGGTCCCCTGGCCCCAGCAGTCATGAGTACCAGGTGCTGCTAGTGTTTCTCTTGCACAGTAATATATGGCCGTGTCCTCGGCTCTCAGGCTGTTCATTTGCAGA
        $total_seq++;
	my $flag = "in-frame";
	my $strand = "Plus";
        my @F=split /\s+/;
	
	# 1.1 . filter the sequence witout V or J gene
	if($F[0] eq "NA" or $F[20] eq "NA")
	{
		$discard{"1"}+=1;
		print E "1\t$ori_line\n";
		$without_vj++;
		next;
	}

	# 1.2. filter the sequence with V and J starand conflict
	if(($F[6]-$F[7])*($F[26]-$F[27])<0)
	{
		$discard{"2"}+=1 and print E "2\t$ori_line\n";
		$vj_conflict++;
		next;
	}

	# 1.3. Convert the sequence and alignment record to plus strand
	my $q_len=length($F[31]);
	my $tmp;
      	if(($F[6]-$F[7])>0){
                $strand="Minus";
                $F[31]=reverse $F[31];
                $F[31]=~tr/ATCGatcg/TAGCtagc/;
                $F[4]=$q_len-$F[4]+1;
                $F[5]=$q_len-$F[5]+1;
                $F[24]=$q_len-$F[24]+1;
                $F[25]=$q_len-$F[25]+1;

                $tmp=$F[4];$F[4]=$F[5];$F[5]=$tmp;
                $tmp=$F[24];$F[24]=$F[25];$F[25]=$tmp;
                $tmp=$F[6];$F[6]=$F[7];$F[7]=$tmp;
                $tmp=$F[26];$F[26]=$F[27];$F[27]=$tmp;
		if($F[10] ne "NA"){
			$F[14]=$q_len-$F[14]+1;
			$F[15]=$q_len-$F[15]+1;
			$tmp=$F[14];$F[14]=$F[15];$F[15]=$tmp;
			$tmp=$F[16];$F[16]=$F[17];$F[17]=$tmp;
		}
        }
 

	# 1.4. Primer Found and correction
	if(defined($inb))
	{
		my $flag_p;
		($flag_p ,@F) = &seq_cut_primer($primer_l_restrict , $p_mis_n ,@F);
		if($flag_p != 2){
			$primer_filter++;
			$discard{"3"}+=1 and print E "3\t$ori_line\n" and next;
		}
	}


	# 2.1. judge the pesudo-gene
        if($F[0] =~ /_unF/ or $F[20] =~ /_unF/){
			$flag = "non-function";
	}

	# 2.2. Find the CDR3 region && traslate it to amino acid
	my $cdr3_ref_start = (split /\./,$F[0])[-2];
	my $cdr3_ref_end = (split /\./,$F[20])[-2];

	my $cdr3_start = $F[5]-($F[7]-$cdr3_ref_start);
	my $cdr3_end = $F[24] + ($cdr3_ref_end-$F[26]);

	my $cdr3length=$cdr3_end-$cdr3_start+1;
	if($cdr3length<=0)# CDR3 length error
	{
		$cdr3_err++;
		$discard{"4"}+=1;print E "4\t$ori_line\n" and next;
	} # alignment error


	my $ORFdna = substr($F[31],$cdr3_start-1,$cdr3length);

	my $aa="NA";
	my $aa_for_clonotypeAA="NA";
	my $V_phase = ($cdr3_start-1)%3;
	my $nt_for_clonotypeAA = $F[31];
	if($cdr3length%3!=0){
		$flag = "out-of-frame(CDR3_length)";
	}
	else{
		$aa=&dna2aa($ORFdna,0);
		$aa_for_clonotypeAA=&dna2aa($nt_for_clonotypeAA,$V_phase);
		$flag = "out-of-frame(stop_codon)" if($aa_for_clonotypeAA=~/\*/);
	}
	
	# store used for sequence abundance fitler
	$cdr3_abund{$ORFdna}++ if($filter_cdr3);
	$clone_abund{$nt_for_clonotypeAA}++ if($filter_clone);


	# store for sequence abundance filter
	
	# 2.3. identify sequence insert/deletion
	# 2.3.1 for Insertion/ D deletion
	my ($VDadd,$DJadd,$VJadd)=("NA","NA","NA");
	my ($v_del,$d_del1,$d_del2,$j_del)=("NA","NA","NA","NA");
	my ($v_len,$j_len)=(length $ref{$F[0]}, length $ref{$F[20]});

	if($ind)# for TRB, IGH
	{
		if($F[10] ne "NA"){
			#Addition of nucleotides between VD genes or DJ genes
			($F[14]-$F[5]>1) and $VDadd=substr($F[31],$F[5],$F[14]-$F[5]-1);
			($F[14]-$F[5]==1) and $VDadd=0;
			($F[24]-$F[15]>1) and $DJadd=substr($F[31],$F[15],$F[24]-$F[15]-1);
			($F[24]-$F[15]==1) and $DJadd=0;
			# for D gene deletion
			my $d_len = length $ref{$F[10]};
			$F[16]>1 and $d_del1=substr($ref{$F[10]},0,$F[16]-1);
			$F[16]==1 and $d_del1=0;
			$F[17]<$d_len and $d_del2=substr($ref{$F[10]},$F[17],$v_len-$F[17]);
			$F[17]==$d_len and $d_del2=0;

		}	

	}else # for TRA, IGL, IGK
	{
		($F[24]-$F[5]>1) and $VJadd=substr($F[31],$F[5],$F[24]-$F[5]-1);
		($F[24]-$F[5]==1) and $VJadd=0;
	}

	# 2.3.2 for V/J Deletion
	$F[7]<$v_len and $v_del=substr($ref{$F[0]},$F[7],$v_len-$F[7]);
	$F[7]==$v_len and $v_del=0;
	$F[26]>1 and $j_del=substr($ref{$F[20]},0,$F[26]-1);
	$F[26]==1 and $j_del=0;

	pop @F;
	my $id = pop @F;
	my $new_line = join "\t" , @F;
	if($filter_cdr3 or $filter_clone){
		print T "$id\t$flag\t$F[0]\t$F[10]\t$F[20]\t$cdr3_start\t$cdr3_end\t$ORFdna\t$aa\t$v_del\t$d_del1\t$d_del2\t$j_del\t$VDadd\t$DJadd\t$VJadd\t$strand\t$nt_for_clonotypeAA\t$aa_for_clonotypeAA\t$new_line\n";
	}else{
		print O "$id\t$flag\t$F[0]\t$F[10]\t$F[20]\t$cdr3_start\t$cdr3_end\t$ORFdna\t$aa\t$v_del\t$d_del1\t$d_del2\t$j_del\t$VDadd\t$DJadd\t$VJadd\t$strand\t$nt_for_clonotypeAA\t$aa_for_clonotypeAA\t$new_line\n";
	}

}
close I;
close T;
#print Dumper(\%cdr3_abund);
#3. sequence abundance filter
my $abund_filter = 0;
if($filter_cdr3)
{
	open I, "$out_dir/$sample_name.structure.tmp.o" or die;
	open L, "|gzip >$out_dir/$sample_name.lowabundance.filter.gz" or die;
	print L "#ID\tfunction\tV_ref\tD_ref\tJ_ref\tCDR3_start\tCDR3_end\tCDR3(dna)\tCDR3(aa)\t3'V_del\t5'D_del\t3'D_del\t5'J_del\tVD_ins\tDJ_ins\tVJ_ins\tstrand\tsequence\tamino_acid\talignment_record\n";
	while(<I>){
		my $cdr3 = (split)[7];
		if($cdr3_abund{$cdr3} >= $filter_cdr3){
			print O "$_";
		}else{
			$abund_filter++;
			print L "$_";
		}
	}
	close I;
	close L;
	unlink "$out_dir/$sample_name.structure.tmp.o";
}
elsif($filter_clone)
{
	open I, "$out_dir/$sample_name.structure.tmp.o" or die;
	open L, "|gzip >$out_dir/$sample_name.lowabundance.filter.gz" or die;
	print L "#ID\tfunction\tV_ref\tD_ref\tJ_ref\tCDR3_start\tCDR3_end\tCDR3(dna)\tCDR3(aa)\t3'V_del\t5'D_del\t3'D_del\t5'J_del\tVD_ins\tDJ_ins\tVJ_ins\tstrand\tsequence\tamino_acid\talignment_record\n";
	while(<I>){
		my $clone = (split)[17];
		if($clone_abund{$clone} >= $filter_clone){
			print O "$_";
		}else{
			$abund_filter++;
			print L "$_";
		}	
	}
	close I;
	unlink "$out_dir/$sample_name.structure.tmp.o";
}
else{
	unlink "$out_dir/$sample_name.structure.tmp.o";
}
close O;


#-----------	out put  ~~~~~~~~~~~~~~~~~~~~~~~~~
print "all_seq: $total_seq\n";
print "Retained_seq: ",$total_seq-$without_vj-$vj_conflict-$cdr3_err-$abund_filter-$primer_filter,"\n";
print "without_V_or_J : $without_vj\n";
print "V_and_J_strand_conflict: $vj_conflict\n";
print "CDR3_less_than_0bp: $cdr3_err\n";
print "seq_abundance_filter: $abund_filter\n";
print "Primer_filter: $primer_filter\n";




# -------------#
# translate 
#--------------#

sub dna2aa{
	my $seq=shift;
#number of nucleotide before condon in the $seq
	my $phase=shift;

	my $aa="";
	$phase==1 and $seq=~s/^\w// and $aa="x";
	$phase==2 and $seq=~s/^\w\w// and $aa="x";
	
	for (my $start=0;$start<length($seq);$start+=3){
		my $codon=uc(substr($seq,$start,3));
		!exists $codon2aa{"$codon"} and $codon2aa{"$codon"}="x";
		$aa.=$codon2aa{"$codon"};
	}
	return $aa;

}


#-----------------------#
#	cut primer seq	#
#-----------------------#
#
# 1. don't cut primer: 1). V/J don't have align result(NA); 
#                      2). V and J align to different strands;
#		       3). error align: forward strand (J at the left of V)/ minus strand(V at the left of J)
# 2. return "NA	NA ...": V or J aligned region all during the primer region

sub seq_cut_primer
{
	my  ($p_l_restrict , $p_mis_n , @line) = @_;
	return(2, @line) if($p_l_restrict==0);
	my ($v_p_len , $j_p_len , $j_len );
	my $primer_flag = 0;
	
	# for V gene

	return (0 , @line) if($line[0] eq "NA" && $line[20] eq "NA"); # seq don't align to VJ ref
	return (0 , @line) if($line[0] ne "NA" && $line[20] ne "NA" && ($line[6]-$line[7])*($line[26]-$line[27])<0);# seq have VJ but align to different strand


	my ($r_start , $r_end , $s_start , $s_end);
	# align to forward strand	***********************************************
	if(($line[0] ne "NA" && $line[6]<$line[7]) or ($line[20] ne "NA" && $line[26]<$line[27]))# align to forward strand
	{
		return (0 , @line) if($line[0] ne "NA" && $line[20] ne "NA" && $line[4]>=$line[24]);# J at the left of V: error and return
	
		# 1 step : conduct V gene
		if($line[0] ne "NA")	
		{
			# caculate the raw postion of primer
			$r_end = (split /\./,$line[0])[3];
			$s_end = $line[4]-($line[6]-$r_end);
			$r_start = 1;
			$s_start = $line[4] - $line[6]+1;

			# primer length restrict
			if($s_end >= $primer_l_restrict && $s_end < length($line[-1]))# 
			{
				# modify position
				if($s_start<=0){
					$s_start = 1;
					$r_start = $line[6]-$line[4]+1;
				}
				# primer check
				my $flag = &primer_align(substr($line[-1],$s_start-1,$s_end-$s_start+1) , substr($ref{$line[0]},$r_start-1,$s_end-$s_start+1) , $p_mis_n);
				 # error correction
				if($flag){
					$primer_flag++;
					substr($line[-1],$s_start-1,$s_end-$s_start+1) = substr($primer_tem_check{$line[0]},$r_start-1,$s_end-$s_start+1);
				}
			}
		}

		# 2 step : conduct J gene
		if($line[20] ne "NA")
		{
			my $raw_s = reverse $primer_tem_check{$line[20]};
			$raw_s =~ tr/ACGT/TGCA/;

			my $p_len = (split /\./,$line[20])[3];
			$r_end =length($ref{$line[20]});
			$r_start = $r_end-$p_len+1;
			$s_start = $line[25]+$r_start-$line[27];
			$s_end = $line[25]+$r_end-$line[27];
			#
			if(length($line[-1])-$s_start+1 >= $primer_l_restrict && $s_start > 0)
			{
				if($s_end>length($line[-1])){
					$s_end = length($line[-1]);
					$r_end = $line[27]+$s_end-$line[25];
				}
				my $flag = &primer_align(substr($line[-1],$s_start-1,$s_end-$s_start+1) , substr($ref{$line[20]},$r_start-1,$s_end-$s_start+1) , $p_mis_n);my $t1=substr($line[-1],$s_start-1,$s_end-$s_start+1);my $t2=substr($ref{$line[20]},$r_start-1,$s_end-$s_start+1);
				if($flag){
					$primer_flag++;
					substr($line[-1],$s_start-1,$s_end-$s_start+1) = substr($raw_s,0,$s_end-$s_start+1);
				}
			}

		}

	}

	

	# align to minus strand		**************************************************
	else # align to minus strand
	{
		return (0 , @line) if($line[0] ne "NA" && $line[20] ne "NA" && $line[4]<=$line[24]);# J at the right of V: error and return
		
		# 1 step : conduct V gene
		if($line[0] ne "NA")
		{	
			my $raw_s = reverse $primer_tem_check{$line[0]};
			$raw_s =~ tr/ACGT/TGCA/;
			my $p_len = length($raw_s);
			
			# caculate the raw postion of primer
			$r_end = 1;
			$r_start = $p_len;
			$s_end = $line[5]+$line[7]-1;
			$s_start = $s_end-$p_len+1;
			
			if($s_start>0 && length($line[-1])-$s_start >= $primer_l_restrict)# primer length filter
			{
				# modify the position
				if($s_end > length($line[-1])){
					$s_end = length($line[-1]);
					$r_end = $line[7]-($s_end-$line[5])
				}

				# primer check
				my $r_sub_s = substr($ref{$line[0]},$r_end-1,$s_end-$s_start+1);
				$r_sub_s = reverse $r_sub_s;
				$r_sub_s =~ tr/ACGT/TGCA/;
				my $flag = &primer_align(substr($line[-1],$s_start-1,$s_end-$s_start+1) , $r_sub_s , $p_mis_n);
				# error correction
				if($flag){
					$primer_flag++;
					substr($line[-1],$s_start-1,$s_end-$s_start+1) = substr($raw_s,0,$s_end-$s_start+1);
				}
				
			}
		}

		# 2 step : conduct J gene
		if($line[20] ne "NA")
		{
			my $raw_s = $primer_tem_check{$line[20]};
			my $p_len = length($raw_s);
			
			# caculate the raw postion of primer
			$r_start = length($ref{$line[20]});
			$r_end = $r_start-$p_len+1;
			$s_start = $line[24]-($r_start-$line[26]);
			$s_end = $line[24]-($r_end-$line[26]);
			
			if($s_end <length($line[-1]) && $s_end >= $primer_l_restrict)# primer length filter
			{
				# modify the position
				if($s_start <= 0){
					$s_start = 1;
					$r_start = $line[26]+$line[24]-1;
				}
				# primer check
				my $r_sub_s = substr($ref{$line[20]},$r_end-1,$s_end-$s_start+1);
				$r_sub_s = reverse $r_sub_s;
				$r_sub_s =~ tr/ACGT/TGCA/;
				my $flag = &primer_align(substr($line[-1],$s_start-1,$s_end-$s_start+1) , $r_sub_s , $p_mis_n);
				# error correction
				if($flag){
					$primer_flag++;
					substr($line[-1],$s_start-1,$s_end-$s_start+1) = substr($raw_s,length($ref{$line[20]})-$r_start,$s_end-$s_start+1);
				}

			}
		}
	}

	return ($primer_flag ,@line);
}

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
