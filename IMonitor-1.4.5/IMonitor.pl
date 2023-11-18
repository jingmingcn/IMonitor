#!/usr/bin/perl -w
use strict;

=head1 Description
        This script use to analyze the immune repertoire sequenced by high throughtput sequencing.
	Input: Pair-end reads with FASTQ format
	       Single-end read with FA format
	       Single-end read with FQ format

=head1 Version
        Version: 1.4.1    Author: wzhang287-c@my.cityu.edu.hk
	Update: 2020.07
		We have added the parameter for the adapter sequence and the discard cutoff
	Version: 1.4.3
	Update: 2021.07
		add the parameters: -addR1/-addR2
	Version: 1.4.4
	Update: 2021.08
		1. update all Ref(Human) and combined chains: TRAB,TRDG,TR,IG
		2. add the function to analyze all chains together
	Version: 1.4.5
	Update: 202109
		1. if CDR3 is too long(there might be two J genes), the J closest to CDR3 will be chose: add the program: Re-alignment.doubleJ.pl

=head1 [parameters]

        -a	<S> full path of input fq file 1
        -b	<S> full path of input fq file 2
	-i	<S> single reads with FASTA format file
	-iq	<S> single reads with FASTQ format file
	-o      <S> output directory path
        -n	<S> sample name, used for prefix of ouput file
	-t	<S> gene type. e.g. TRA,TRB,TRD,TRG,IGH,IGK,IGL,IGKL,TRAB,TRDG,TR,IG,All. here, TRAB=TRA+TRB, TRDG=TRD+TRG,TR=TRA+TRB+TRD+TRG,IG=IGH+IGKL 
	-k      <I> read length [100]

	-d	<F> add the paremeter means consider D genes for analysis. For IGH,TRB is necessary
	-c	    logical value, for without V or J sequence, find the CDR3 by conservative region. V(YXC),J([WF]GXG)
	-jif	<I> J alignment identity for filtering [80]
	-vif	<I> V alignment identity for filtering [80]
	-r      <S> The reference directory [Bin/Ref/gene_type]
	-Q	<I> sequencing quality for filtering [15]
	-RQ	<Fl> qulity filter rate, used with -Q [0.05]
	-f1	<I> CDR3(nucleotide) abundance filter [0]
	-f2     <I> Clonotype(nucleotide, ful-length) abundance filter [0]
	-m	    logical value, used to analyze hyper-mutation
	-ew	    logical value, sequencing error correction for whole sequence. but this need a long time to run,only for FASTQ files as input
	-ec	    logical value, sequencing error correction for only CDR3.only for FASTQ files as input
	-Qe	<I> the quality(as a cutoff) is used for sequencing error correction. [20](-ew or -ec is required) 
	-adseq	<S>	3' adapter sequence
	-adcut	<I>	cut sequence from adaptor index,unless performed -f/-r also in use
				discard the read when the adaptor index of the read is less than INT
	-addR1	    logical value, add read1 for further analysis if the PE reads are failed to merge
	-addR2	    logical value, add read2 for further analysis if the PE reads are failed to merge

	the next two options only for sequencing type
	-v	<I> used to calculate the base quality [64]
	-seqType <I> Sequence fq type, 0->old fastq(Hiseq2000), 1->new fastq(Hiseq4000)[default: 0]
	old fastq id: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2
	new fastq id: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC
	
	-mul	<I> split the fa into multiple parts for alignment [3]
	-Rs	the R script directory.for a new system, this parameter is need to change[/ifs1/ST_MED/USER/zhangwei/software/R-3.0.2/bin/Rscript]
	-h	print help information
	

=head1 Note

	#1. If  Pair-end(PE) sequencing FASTQ format as input, then:
		perl IMonitor.pl
		Compulsory: -a -b -o -n -t -Rs
		Optionally: others
	#2. If Single-end(SE) sequencing FASTA format as input, then:
		perl IMonitor.pl
		Compulsory: -i -o -n -t -Rs
		Optionally: others, but -ew,-ec are invalid here
	#3. If Single-end(SE) sequencing FASTQ format as input, then:
		perl IMonitor.pl
		Compulsory: -iq -o -n -t -Rs
		Optionally: others
	    For Zebra sequencing(Single-end), Compulsory: -iq -o -n -t -Rs -v 33 [-Qe 25 (-ew or -ec is required)]
	#4. If Pair-end(PE) sequencing FASTQ format as input, and add the single reads(the PE reads are failed to merge), then:
		perl IMonitor.pl
		Compulsory: -a -b -o -n -t -addR1(or -addR2)
	#5. for the raw sequencing data including multiple chains, the paramter -t (TRAB,TRDG,TR,IG) can be used for analysis together
		
	# Note:
	  The rate of IMonitor output is multipled 100%


=cut
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);

my ($fq_1,$fq_2,$ad1, $ad2, $get_left_fq,$refdir,$fa_ref,$outdir,$basename,$gene_type,$primer_f,$qsub,$D_prob,$readLength,$filter_Q,$filter_rate,$filter_cdr3,$filter_clone,$mutation_f,$error_f,$satruation_f,$Help,$Rscript,$cal_qual_v,$jif,$vif,$quene,$mul_num,$adseq,$adcut);
$D_prob = 0 ; # default
my ($error_f_c , $single_file, $cdr3_byconserve, $single_file_q, $qual_error_c, $add_read1_f, $add_read2_f);
$cdr3_byconserve = 0;
($add_read1_f, $add_read2_f) = (0,0);
my $seqtype;

GetOptions(
        "a:s" => \$fq_1,
        "b:s" => \$fq_2,
	#	"A1:s" => \$ad1,
	#	"A2:s" => \$ad2,
	"i:s" =>\$single_file,
	"iq:s" => \$single_file_q,
        "r:s" => \$refdir,
        "o:s" => \$outdir,
        "n:s" => \$basename,
	"t:s" => \$gene_type,
	"d" => \$D_prob,
	"c" => \$cdr3_byconserve,
	"jif:i" => \$jif,
	"vif:i" => \$vif,
	"k:i" => \$readLength,
	"Q:i" => \$filter_Q,
	"RQ:f" => \$filter_rate,
	"f1:i" => \$filter_cdr3,
	"f2:i" => \$filter_clone,
	"m" => \$mutation_f,
	"ew" => \$error_f,
	"ec" => \$error_f_c,
	"Qe" => \$qual_error_c,
	"adseq:s" => \$adseq,
	"adcut:i" => \$adcut,
	"mul:i" => \$mul_num,
	"v:i" => \$cal_qual_v,
	"Rs:s" => \$Rscript,
	"seqType:i" => \$seqtype,
	"addR1" => \$add_read1_f,
       	"addR2" => \$add_read2_f,
        "help!" => \$Help,

);

die `pod2text $0` if ($Help or (!$single_file_q and !$single_file and (!$fq_1 or !$fq_2)) or !$gene_type or !$outdir or !$basename);


#-------- give the dedault parameter	

$readLength||=100;
$adseq||="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA";
$adcut||=80;
$filter_Q=15 unless(defined($filter_Q));
$filter_rate=0.05 unless(defined($filter_rate));
$filter_cdr3=0 unless(defined($filter_cdr3));
$filter_clone=0 unless(defined($filter_clone));
$cal_qual_v=64 unless(defined($cal_qual_v));
$jif=80 unless(defined($jif));
$vif=80 unless(defined($vif));
$mul_num=3 unless(defined($mul_num));
$satruation_f = 1;
$seqtype = 0 unless(defined($seqtype));
$qual_error_c = 20 unless(defined($qual_error_c));

$Rscript = "/data/Public_tools/R-4.0.2/bin/Rscript" unless(defined($Rscript));



#----- change to absolute path
my $current_path=$ENV{"PWD"};

my ($fq1_base, $fq2_base, $ad1_dir);
if(defined($fq_1)){
	$fq1_base=basename $fq_1;
	$fq2_base=basename $fq_2;
	$ad1_dir=dirname $fq_1;
	$fq_1="$current_path/$fq_1" unless($fq_1=~/^\//);
	$fq_2="$current_path/$fq_2" unless($fq_2=~/^\//);
	#	$ad1="$current_path/$ad1" unless($ad1=~/^\//);
	#	$ad2="$current_path/$ad2" unless($ad2=~/^\//);
}elsif(defined($single_file)){
	$single_file = "$current_path/$single_file" unless($single_file=~/^\//)
}elsif(defined($single_file_q)){
	$single_file_q = "$current_path/$single_file_q" unless($single_file_q=~/^\//)
}
$outdir="$current_path/$outdir" unless($outdir=~/^\//);
$refdir="$current_path/$refdir" if(defined($refdir) and !($refdir=~/^\//));



#------- the default reference 
if(!defined($refdir) || $refdir eq "")
{
	if(($gene_type eq "IGH" || $gene_type eq "IGKL" || $gene_type eq "TRA" || $gene_type eq "TRB" || $gene_type eq "TRD" || $gene_type eq "TRG" || $gene_type eq "IGK" || $gene_type eq "IGL" || $gene_type eq "TRAB" || $gene_type eq "TRDG" || $gene_type eq "TR" || $gene_type eq "IG" || $gene_type eq "All"))
	{
		$refdir = "$Bin/Ref/$gene_type";
	}else
	{
		print "reference dir error!\n";
	}
}


#--- make some directory
system "mkdir -p $outdir/Figures" unless(-d "$outdir/Figures");
system "mkdir -p $outdir/Result" unless(-d "$outdir/Result");
system "mkdir -p $outdir/Bin" unless(-d "$outdir/Bin");
if(-d "$outdir/Align"){unlink glob "$outdir/Align/*.sh";
}else{system "mkdir -p $outdir/Align";}


my $u_argu_for_cmr=$readLength-1;
my $min_FqMergerJava=$u_argu_for_cmr;
my $max_FqMergerJava=$readLength*2-50;


open EX, ">$outdir/Bin/$basename.Execute_all.sh" or die;

# 1. raw fastq sequence process	---------------

open O, ">$outdir/Bin/$basename.merge_fq_fq2fa.sh" or die;
print O "#!/bin/sh\n\n";
print O "echo \"$basename.merge_fq_fq2fa.sh start: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
if(defined($single_file))
{
	print O "# split fa sequence ----------\n";
	print O "perl $Bin/Merge_filter_fa.pl $single_file $outdir/Result/$basename.merged.fa.gz $outdir/Result/$basename.change_id.backup.gz $outdir/Result/${basename}_insert_size_len > $outdir/${basename}_merge_filter_1.txt\n";
	print O "perl $Bin/Split_fq_into_part.pl $outdir/Result/$basename.merged.fa.gz $outdir/${basename}_merge_filter_1.txt $outdir/Align $basename $mul_num\n";
}elsif(defined($single_file_q)){
	print O "# low quality and split fa sequence ----------\n";
	print O "$Bin/SOAPnuke filter -l 20 -Q 2 -5 0 -q 0.1 -m 25 -n 0.05 -p 1 -G -1 $single_file_q -o $outdir -f $adseq --cutAdaptor $adcut -C $outdir/${basename}_filter_1.fq.gz\n";
	print O "perl $Bin/change_format.pl $outdir/Basic_Statistics_of_Sequencing_Quality.txt $outdir/${basename}_raw_data_filter_1.txt\n";
	print O "perl $Bin/Single_read_FQ_tofa.pl $outdir/${basename}_filter_1.fq.gz $filter_Q $filter_rate $outdir/Result/$basename.merged.fa.gz $outdir/Result/$basename.change_id.backup.gz $outdir/Result/${basename}_insert_size_len $cal_qual_v >$outdir/${basename}_merge_filter_1.txt\n\n";
	print O "perl $Bin/Split_fq_into_part.pl $outdir/Result/$basename.merged.fa.gz $outdir/${basename}_merge_filter_1.txt $outdir/Align $basename $mul_num\n";
}else{

	print O "# 1. low quality and adapter filter--------------\n";
	print O "perl $Bin/Filter_adapter_lowqual_N.pl $fq_1 $fq_2 $readLength $outdir $basename $seqtype -v $cal_qual_v >$outdir/${basename}_raw_data_filter_1.txt\n";
	print O "# 2. merge PE read-------------\n";
	print O "$Bin/src/cope -a $outdir/${basename}_filter_1.fq.gz -b $outdir/${basename}_filter_2.fq.gz -o $outdir/$basename.merged_fq -2 $outdir/$fq1_base.left -3 $outdir/$fq2_base.left -l 10 -u $u_argu_for_cmr -c 0.9 -m 0 -s $cal_qual_v >$outdir/$basename.fqMerging.log 2>$outdir/$basename.fqMerging.err\n";
	print O "java -jar -Xmx2G $Bin/FqMerger.jar $outdir/$fq1_base.left $outdir/$fq2_base.left $outdir/$basename.further_merged.gz 50 $max_FqMergerJava $readLength 0.9 0.7\n";

	print O "# 3. filter and change into fa--------------\n";
	if($add_read1_f or $add_read2_f)# add the single read if PE reads failed to merge
	{
		my $single_f_n = "$outdir/$fq1_base.left";
		$single_f_n = "$outdir/$fq2_base.left" if($add_read2_f);
		print O "perl $Bin/Merge_filter_fqtofa.V2.pl $outdir/$basename.merged_fq $outdir/$basename.further_merged.gz $outdir/$basename.fqMerging.log $filter_Q $filter_rate $outdir/Result/$basename.merged.fa.gz $outdir/Result/$basename.change_id.backup.gz $outdir/Result/${basename}_insert_size_len $cal_qual_v $single_f_n >$outdir/${basename}_merge_filter_1.txt\n";
		print O "mv $outdir/$basename.further_merged.gz.gz $outdir/$basename.further_merged.gz\n\n";
	}else{
		print O "perl $Bin/Merge_filter_fqtofa.V2.pl $outdir/$basename.merged_fq $outdir/$basename.further_merged.gz $outdir/$basename.fqMerging.log $filter_Q $filter_rate $outdir/Result/$basename.merged.fa.gz $outdir/Result/$basename.change_id.backup.gz $outdir/Result/${basename}_insert_size_len $cal_qual_v >$outdir/${basename}_merge_filter_1.txt\n\n";
	}
	print O "# 4. split fa sequence ----------\n";
	
	print O "perl $Bin/Split_fq_into_part.pl $outdir/Result/$basename.merged.fa.gz $outdir/${basename}_merge_filter_1.txt $outdir/Align $basename $mul_num\n";

}

print O "echo \"$basename.merge_fq_fq2fa.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
close O;
print EX "sh $outdir/Bin/$basename.merge_fq_fq2fa.sh\n";

# 2. Blast aligment && Re-alingment	----------------------
my %jobid;
my @last_job_id_2;
#We use different -W arguements for different types of genes
my %W_argu=("V"=>15, "D"=>4, "J"=>10);
my %v_c=("V"=>1, "D"=>3, "J"=>1);
my %b_c=("V"=>3, "D"=>5, "J"=>3);

for(my $i=1 ; $i<=$mul_num ; $i++)
{
	open O, ">$outdir/Bin/$basename.blast.$i.sh";
	print O "#!/bin/sh\n\n";
	print O "echo \"$basename.blast.$i.sh start: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
	for my $gene('V','D','J')
	{
		next if(!$D_prob and $gene eq "D");
	
		print O "$Bin/blastall -p blastn -d $refdir/$gene_type.$gene.ref_index -i $outdir/Align/$basename.fa.tmp.$i -o $outdir/Align/$basename.blast.$gene.m8.$i -y 15 -r 5 -q -4 -g F -W $W_argu{$gene} -K 3 -v $v_c{$gene} -b $b_c{$gene} -m 8 -a 6\ngzip -f $outdir/Align/$basename.blast.$gene.m8.$i\n";
		
	}
	if($D_prob){
		print O "perl $Bin/Re-alignment.pl $outdir/Align/$basename.fa.tmp.$i $outdir/Align/$basename.blast.V.m8.$i.gz $outdir/Align/$basename.blast.J.m8.$i.gz  $outdir/Align/$basename.blast.ref.SingleLine.$i $basename $refdir/$gene_type.converted.fa $gene_type -F $outdir/Align/$basename.blast.D.m8.$i.gz -jif $jif -vif $vif >$outdir/Align/${basename}_VDJ.alignment.stat_2.$i.txt\n";
		print O "perl $Bin/Re-alignment.doubleJ.pl $outdir/Align/$basename.blast.ref.SingleLine.$i $outdir/Align/$basename.blast.V.m8.$i.gz $outdir/Align/$basename.blast.J.m8.$i.gz  $outdir/Align/$basename.blast.ref.SingleLine.$i.tmp $basename $refdir/$gene_type.converted.fa $gene_type -F $outdir/Align/$basename.blast.D.m8.$i.gz -jif $jif -vif $vif >$outdir/Align/${basename}_VDJ.alignment.stat_2.$i.txt\n";
		print O "mv $outdir/Align/$basename.blast.ref.SingleLine.$i.tmp $outdir/Align/$basename.blast.ref.SingleLine.$i\n";
	}else{
		print O "perl $Bin/Re-alignment.pl $outdir/Align/$basename.fa.tmp.$i $outdir/Align/$basename.blast.V.m8.$i.gz $outdir/Align/$basename.blast.J.m8.$i.gz  $outdir/Align/$basename.blast.ref.SingleLine.$i $basename $refdir/$gene_type.converted.fa $gene_type  -jif $jif -vif $vif >$outdir/Align/${basename}_VDJ.alignment.stat_2.$i.txt\n";
		print O "perl $Bin/Re-alignment.doubleJ.pl $outdir/Align/$basename.blast.ref.SingleLine.$i $outdir/Align/$basename.blast.V.m8.$i.gz $outdir/Align/$basename.blast.J.m8.$i.gz  $outdir/Align/$basename.blast.ref.SingleLine.$i.tmp $basename $refdir/$gene_type.converted.fa $gene_type  -jif $jif -vif $vif >$outdir/Align/${basename}_VDJ.alignment.stat_2.$i.txt\n";
		print O "mv $outdir/Align/$basename.blast.ref.SingleLine.$i.tmp $outdir/Align/$basename.blast.ref.SingleLine.$i\n";
	}
	print O "echo \"$basename.blast.$i.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
	close O;
	print EX "sh $outdir/Bin/$basename.blast.$i.sh\n";

}

# 3. re-alignment and structure identify	--------------------------


open O, ">$outdir/Bin/$basename.structure.sh" or die;
print O "#!/bin/sh\n\n";
print O "echo \"$basename.structure.sh start: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
#--------------- cat the aligment result --------------------
print O "# 1. cat the V(D)J alignment and re-alignment------------\n";

my $cat_singleline = "cat $outdir/Align/$basename.blast.ref.SingleLine.1";
for(2..$mul_num){$cat_singleline= "$cat_singleline $outdir/Align/$basename.blast.ref.SingleLine.$_"}

print O "$cat_singleline |gzip -cf>$outdir/$basename.blast.ref.SingleLine.gz\n";
print O "perl $Bin/Cat_alignment_stat.pl $outdir/Align $basename $outdir/${basename}_VDJ.alignment.stat_2.txt $mul_num\n";

#unlink("$outdir/Align/qsub_sh_all.e") if(-e "$outdir/Align/qsub_sh_all.e");
#unlink("$outdir/Align/qsub_sh_all.o") if(-e "$outdir/Align/qsub_sh_all.o");
#print O "head -20 $outdir/Align/*.sh.e* >$outdir/Align/qsub_sh_all.e\nhead -20 $outdir/Align/*.sh.o* >$outdir/Align/qsub_sh_all.o\n\n";
#print O "rm $outdir/Align/$basename.blast.ref.SingleLine.[1-$mul_num] $outdir/Align/${basename}_VDJ.alignment.stat_2.[1-$mul_num].txt $outdir/Align/*.sh.e* $outdir/Align/*.sh.o*\n\n";

# PCR and sequencing error correction
my $singleline = "$outdir/$basename.blast.ref.SingleLine.gz";


print O "perl $Bin/Change_combine_seq_id.pl $outdir/$basename.blast.ref.SingleLine.gz $outdir/Result/$basename.change_id.backup.gz $outdir/$basename.blast.ref.SingleLine.new.gz\n";
print O "mv $outdir/$basename.blast.ref.SingleLine.new.gz $outdir/$basename.blast.ref.SingleLine.gz\n";

my $fq_all = "$outdir/$basename.merged_fq.all";
if(defined($single_file_q)){
	$fq_all = "$outdir/${basename}_filter_1.fq.gz";
}

if($error_f)
{
	print O "# 2. PCR and sequencing error correction----------\n";
	unless($single_file_q){
                if($add_read1_f or $add_read2_f){
                        my $single_f_n = "$outdir/$fq1_base.left";
                        $single_f_n = "$outdir/$fq2_base.left" if($add_read2_f);
			print O "gzip -dc $outdir/$basename.further_merged.gz>$outdir/$basename.further_merged\ncat $outdir/$basename.merged_fq $outdir/$basename.further_merged $single_f_n > $fq_all\n";
		}else{
			print O "gzip -dc $outdir/$basename.further_merged.gz>$outdir/$basename.further_merged\ncat $outdir/$basename.merged_fq $outdir/$basename.further_merged > $fq_all\n";
		}
	}
	print O "perl $Bin/Filter_for_correct.pl $refdir/$gene_type.converted.fa $outdir/$basename.blast.ref.SingleLine.gz $fq_all $outdir/Result/$basename.change_id.backup.gz $outdir/$basename.effect.data.gz >$outdir/${basename}_pcr_seq_error_cor_3.txt\n";
	print O "perl $Bin/Correct_seq_pcr_error.pl -i $outdir/$basename.effect.data.gz -d $outdir/Result/${basename}_seq_error_discard.gz -o $outdir/$basename.error.correct.gz -M2 0 -v $cal_qual_v -Q $qual_error_c >>$outdir/${basename}_pcr_seq_error_cor_3.txt\n";
	$singleline = "$outdir/$basename.error.correct.gz";
}
elsif($error_f_c)# Correct only for CDR3 region
{
	print O "# 2. PCR and sequencing error correction----------\n";
	unless($single_file_q){
                if($add_read1_f or $add_read2_f){
                        my $single_f_n = "$outdir/$fq1_base.left";
                        $single_f_n = "$outdir/$fq2_base.left" if($add_read2_f);
			print O "gzip -dc $outdir/$basename.further_merged.gz>$outdir/$basename.further_merged\ncat $outdir/$basename.merged_fq $outdir/$basename.further_merged $single_f_n > $fq_all\n";
		}else{
			print O "gzip -dc $outdir/$basename.further_merged.gz>$outdir/$basename.further_merged\ncat $outdir/$basename.merged_fq $outdir/$basename.further_merged > $fq_all\n";
		}
	}
	print O "perl $Bin/Filter_for_correct_CDR3.pl $refdir/$gene_type.converted.fa $outdir/$basename.blast.ref.SingleLine.gz $fq_all $outdir/Result/$basename.change_id.backup.gz $outdir/$basename.effect.data.gz >$outdir/${basename}_pcr_seq_error_cor_3.txt\n";
	print O "perl $Bin/Correct_seq_pcr_error.pl -i $outdir/$basename.effect.data.gz -d $outdir/Result/${basename}_seq_error_discard.gz -o $outdir/$basename.error.correct.gz -c -M1 3 -M2 0 -v $cal_qual_v -Q $qual_error_c >>$outdir/${basename}_pcr_seq_error_cor_3.txt\n";
	$singleline = "$outdir/$basename.error.correct.gz";
}

print O "# 3. identify the sequence structure----------------\n";
if($D_prob){
	print O "perl $Bin/Structural_analysis.pl -i $singleline -o $outdir/Result -n $basename -r $refdir/$gene_type.converted.fa -d -f1 $filter_cdr3 -f2 $filter_clone >$outdir/${basename}_structure.stat_3.txt\n";
}else{
	print O "perl $Bin/Structural_analysis.pl -i $singleline -o $outdir/Result -n $basename -r $refdir/$gene_type.converted.fa -f1 $filter_cdr3 -f2 $filter_clone >$outdir/${basename}_structure.stat_3.txt\n";
}

#----   find CDR3 by conserved amino acid

if($cdr3_byconserve && $error_f_c)
{
        print O "perl $Bin/Find_CDR3_by_conserve.pl $outdir/$basename.blast.ref.SingleLine.gz $outdir/$basename.structure.cdr3.by.conserve  >/dev/null\n";
	print O "perl $Bin/Filter_for_correct_CDR3.conserve.pl $refdir/$gene_type.converted.fa $outdir/$basename.structure.cdr3.by.conserve $fq_all $outdir/Result/$basename.change_id.backup.gz > $outdir/$basename.filter.conserve\n";
	print O "perl $Bin/Correct_seq_pcr_error.pl -i $outdir/$basename.filter.conserve -d $outdir/Result/${basename}_seq_error_discard.conserve.gz -o $outdir/$basename.error.correct.conserve.gz -c -M1 3 -M2 0 -v $cal_qual_v -Q $qual_error_c >>$outdir/${basename}_pcr_seq_error_cor_3.txt\n";	
	print O "perl $Bin/Cat_error_stat.pl $outdir/${basename}_pcr_seq_error_cor_3.txt >$outdir/o\nmv $outdir/o $outdir/${basename}_pcr_seq_error_cor_3.txt\n";
	print O "perl $Bin/Find_CDR3_by_conserve.pl $outdir/$basename.error.correct.conserve.gz $outdir/$basename.structure.cdr3.by.conserve  >>$outdir/${basename}_structure.stat_3.txt\n";
        print O "mv $outdir/Result/$basename.structure.gz $outdir/$basename.structure.1.gz\nzcat $outdir/$basename.structure.1.gz|cat - $outdir/$basename.structure.cdr3.by.conserve >$outdir/Result/$basename.structure\ngzip -f $outdir/Result/$basename.structure\n";
}
elsif($cdr3_byconserve && $error_f)
{
	print O "perl $Bin/Find_CDR3_by_conserve.pl $outdir/$basename.blast.ref.SingleLine.gz $outdir/$basename.structure.cdr3.by.conserve  >/dev/null\n";
	print O "perl $Bin/Filter_for_correct.conserve.pl $refdir/$gene_type.converted.fa $outdir/$basename.structure.cdr3.by.conserve $fq_all $outdir/Result/$basename.change_id.backup.gz > $outdir/$basename.filter.conserve\n";
	print O "perl $Bin/Correct_seq_pcr_error.pl -i $outdir/$basename.filter.conserve -d $outdir/Result/${basename}_seq_error_discard.conserve.gz -o $outdir/$basename.error.correct.conserve.gz -M2 0 -v $cal_qual_v -Q $qual_error_c >>$outdir/${basename}_pcr_seq_error_cor_3.txt\n";
	print O "perl $Bin/Cat_error_stat.pl $outdir/${basename}_pcr_seq_error_cor_3.txt >$outdir/o\nmv $outdir/o $outdir/${basename}_pcr_seq_error_cor_3.txt\n";
	print O "perl $Bin/Find_CDR3_by_conserve.pl $outdir/$basename.error.correct.conserve.gz $outdir/$basename.structure.cdr3.by.conserve  >>$outdir/${basename}_structure.stat_3.txt\n";
	print O "mv $outdir/Result/$basename.structure.gz $outdir/$basename.structure.1.gz\nzcat $outdir/$basename.structure.1.gz|cat - $outdir/$basename.structure.cdr3.by.conserve >$outdir/Result/$basename.structure\ngzip -f $outdir/Result/$basename.structure\n";
}
elsif($cdr3_byconserve)
{
	print O "perl $Bin/Find_CDR3_by_conserve.pl $outdir/$basename.blast.ref.SingleLine.gz $outdir/$basename.structure.cdr3.by.conserve  >>$outdir/${basename}_structure.stat_3.txt\n";
	print O "mv $outdir/Result/$basename.structure.gz $outdir/$basename.structure.1.gz\nzcat $outdir/$basename.structure.1.gz|cat - $outdir/$basename.structure.cdr3.by.conserve >$outdir/Result/$basename.structure\ngzip -f $outdir/Result/$basename.structure\n";
}


print O "\n";

# for old version of reference, the id need to changed
if(-e "$refdir/$gene_type.backup" and -s "$refdir/$gene_type.backup")
{
	print O "perl $Bin/Change_mpcr_id.pl $outdir/Result/$basename.structure.gz $refdir/$gene_type.backup $outdir/$basename.structure.tmp.gz\n";
	print O "mv $outdir/$basename.structure.tmp.gz $outdir/Result/$basename.structure.gz\n";
}

print O "echo \"$basename.structure.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
close O;
print EX "sh $outdir/Bin/$basename.structure.sh\n";



# 4. statics and draw the figures	--------------------
open O, ">$outdir/Bin/$basename.statistics.graph.sh";
print O "#!/bin/sh\n\n";
print O "echo \"$basename.statistics.graph.sh start: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
my $flag = 1;
$flag =0 if(defined($single_file));
my $error_flag = 0;
$error_flag = 1 if($error_f || $error_f_c);
print O "\nperl $Bin/Get_summary_stat.pl $outdir $basename $flag $error_flag $cdr3_byconserve > $outdir/Result/${basename}_bascial_filter_stat.txt\n\n";

if($gene_type eq "TRAB" || $gene_type eq "TRDG" || $gene_type eq "TR" || $gene_type eq "IG" || $gene_type eq "All")
{
	my $out_sh = "perl $Bin/Statistics.pl -i $outdir/Result/$basename.structure.gz -o $outdir/Result -c $refdir/$gene_type.alignment -r $refdir/$gene_type.backup";
	$out_sh = "$out_sh -d " if($D_prob);
	$out_sh = "$out_sh -m " if($mutation_f);
	$out_sh = "$out_sh -s " if($satruation_f);
	#	$out_sh = "$out_sh >$outdir/Result/${basename}_further.stat.txt";
	
	my $draw_sh = "perl $Bin/R_draw_figure.pl -i $outdir/Result -o $outdir/Figures -f $outdir/Result/${basename}_insert_size_len -c -R $Rscript";
	$draw_sh .= " -s" if($satruation_f);
	$draw_sh .= " -m" if($mutation_f);
	if($gene_type eq "TRAB"){
		print O "$out_sh -t TRA -n ${basename}_TRA >$outdir/Result/${basename}_TRA_further.stat.txt\n";
		print O "$draw_sh -n ${basename}_TRA\n";
		print O "$out_sh -t TRB -n ${basename}_TRB >$outdir/Result/${basename}_TRB_further.stat.txt\n";
		print O "$draw_sh -n ${basename}_TRB\n";
	}elsif($gene_type eq "TRDG"){
                print O "$out_sh -t TRD -n ${basename}_TRD >$outdir/Result/${basename}_TRD_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_TRD\n";
                print O "$out_sh -t TRG -n ${basename}_TRG >$outdir/Result/${basename}_TRG_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_TRG\n";		
	}elsif($gene_type eq "TR"){
                print O "$out_sh -t TRA -n ${basename}_TRA >$outdir/Result/${basename}_TRA_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_TRA\n";
                print O "$out_sh -t TRB -n ${basename}_TRB >$outdir/Result/${basename}_TRB_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_TRB\n";
                print O "$out_sh -t TRD -n ${basename}_TRD >$outdir/Result/${basename}_TRD_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_TRD\n";
                print O "$out_sh -t TRG -n ${basename}_TRG >$outdir/Result/${basename}_TRG_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_TRG\n";
        }elsif($gene_type eq "IG"){
                print O "$out_sh -t IGH -n ${basename}_IGH >$outdir/Result/${basename}_IGH_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_IGH\n";
                print O "$out_sh -t IGKL -n ${basename}_IGKL >$outdir/Result/${basename}_IGKL_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_IGKL\n";
	}elsif($gene_type eq "All"){
                print O "$out_sh -t TRA -n ${basename}_TRA >$outdir/Result/${basename}_TRA_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_TRA\n";
                print O "$out_sh -t TRB -n ${basename}_TRB >$outdir/Result/${basename}_TRB_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_TRB\n";
                print O "$out_sh -t TRD -n ${basename}_TRD >$outdir/Result/${basename}_TRD_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_TRD\n";
                print O "$out_sh -t TRG -n ${basename}_TRG >$outdir/Result/${basename}_TRG_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_TRG\n";
                print O "$out_sh -t IGH -n ${basename}_IGH >$outdir/Result/${basename}_IGH_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_IGH\n";
                print O "$out_sh -t IGKL -n ${basename}_IGKL >$outdir/Result/${basename}_IGKL_further.stat.txt\n";
                print O "$draw_sh -n ${basename}_IGKL\n";
	}

}else{
my $out_sh = "perl $Bin/Statistics.pl -i $outdir/Result/$basename.structure.gz -o $outdir/Result -n $basename -c $refdir/$gene_type.alignment -r $refdir/$gene_type.backup";

$out_sh = "$out_sh -d " if($D_prob);
$out_sh = "$out_sh -m " if($mutation_f);
$out_sh = "$out_sh -s " if($satruation_f);
$out_sh = "$out_sh >$outdir/Result/${basename}_further.stat.txt";

print O "$out_sh\n";

#my $flag = 1;
#$flag =0 if(defined($single_file));
#my $error_flag = 0;
#$error_flag = 1 if($error_f || $error_f_c);
#print O "\nperl $Bin/Get_summary_stat.pl $outdir $basename $flag $error_flag $cdr3_byconserve > $outdir/Result/${basename}_bascial_filter_stat.txt\n\n";

my $draw_sh = "perl $Bin/R_draw_figure.pl -n $basename -i $outdir/Result -o $outdir/Figures -f $outdir/Result/${basename}_insert_size_len -c -R $Rscript";
$draw_sh .= " -s" if($satruation_f);
$draw_sh .= " -m" if($mutation_f);
print O "$draw_sh\n";
}

print O "echo \"$basename.statistics.graph.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
close O;
print EX "sh $outdir/Bin/$basename.statistics.graph.sh\n";

#-------------	remove the process files`--------------
open O, ">$outdir/Bin/$basename.rm.intermediate.file.sh" or die;
print O "#!/bin/sh\n\n";
print O "echo \"$basename.rm.intermediate.file.sh start: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
if(!defined($single_file) and !defined($single_file_q)){
	print O "rm $outdir/${basename}_filter_1.fq.gz $outdir/${basename}_filter_2.fq.gz $outdir/$fq1_base.left $outdir/$fq2_base.left\n";
	print O "rm $outdir/$basename.further_merged* $outdir/$basename.merged_fq\n";
	print O "rm $outdir/$basename.fqMerging.log $outdir/$basename.fqMerging.err $outdir/${basename}_raw_data_filter_1.txt\n";
}
if(defined($single_file_q)){
	print O "rm $outdir/${basename}_raw_data_filter_1.txt $outdir/Basic_Statistics_of_Sequencing_Quality.txt $outdir/Statistics_of_Filtered_Reads.txt $outdir/Base_distributions_by_read_position_1.txt $outdir/Distribution_of_Q20_Q30_bases_by_read_position_1.txt $outdir/Base_quality_value_distribution_by_read_position_1.txt $outdir/${basename}_filter_1.fq.gz\n";
}
print O "rm $outdir/${basename}_merge_filter_1.txt\n";

print O "rm $outdir/Align/*\n";

print O "rm $outdir/${basename}_VDJ.alignment.stat_2.txt $outdir/$basename.blast.ref.SingleLine.gz $outdir/${basename}_structure.stat_3.txt\n";

if($error_f_c || $error_f){
	print O "rm $outdir/$basename.effect.data.gz $outdir/${basename}_pcr_seq_error_cor_3.txt $outdir/$basename.error.correct.gz\n";
	print O "rm $outdir/$basename.merged_fq.all\n" if(!defined($single_file) and !defined($single_file_q));
	print O "rm $outdir/$basename.error.correct.conserve.gz  $outdir/$basename.filter.conserve\n" if($cdr3_byconserve);
}
if($cdr3_byconserve){
	print O "rm $outdir/$basename.structure.cdr3.by.conserve $outdir/$basename.structure.1.gz\n";
}


print O "echo \"$basename.rm.intermediate.file.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
close O;
print EX "sh $outdir/Bin/$basename.rm.intermediate.file.sh\n";

close EX;
