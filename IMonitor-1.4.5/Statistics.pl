#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;


=head1	Funtion
	
	perl Statistics.pl -i -o -n []	

=head1 Parameter
	
	-i	input file
	-o	output directory
	-n	sample name
	-r      reference sequence id file
	-d	logical value, need to analyze the D gene
	-w	a abundance window for count the CDR3(aa) frequency distribution [1000]
	-c	analyze the nucleotide compositon, input the reference with multiple alignment
	-t	chain, eg. TRA, TRB,TRD,TRG,IGH,IGKL,IGK,IGL

	-m	logical value, analyze the hyper-mutation
	-mf	the minimum CDR3 AA abundance which used for hyper-mutation [10]
	

	-s 	logical value, for saturation curve(rarefraction curve)
	-1	the minimum CDR3 AA abundance, used for saturation curve, with -s,[3]

	Set filtered conditions for statistics:
	-nf	whether need to filter non-functional sequences, including pesdo-/stop codon/incorrect CDR3 length
	-cf	whether need to filter sequence that don't contain amino acid "C","[F/W]"
	-lf	whether need to filter very short CDR3 less than 3 (<3)


=head1	Description
	1. sequence mapped to "non-function" is filtered
	2. output V/J/V-J usage
	3. output CDR3/Clone(nt and aa) frequency distribution(files)
	4. output CDR3(nt) length distribution/frequency distribution
	5. output insertion/deletion length distribution
	6. output V/D/J length dis in CDR3 region
	7. output a basical statitics file
	7. optionally, V/J nuleotide composition
	8. optionally, saturation curve(rarefraction curve) data

	201607 add:
	stat CDR3 found by conserative region, without V or J. However, only CDR3-length/CDR3_NT-freq/CDR3_AA_section/Clone_NT-freq add this, others will also statistic sequence with both  V and J

	201710 add:
	-nf -cf -lf some filtered conditions will be used to filter some sequences

=cut

my ($file,$out_dir,$d_flag, $name,$sature_f,$mut_f, $mini_abund_mut,$ref_s , $one , $window, $ref_id);
my ($non_func_f, $cdr3_len_f, $aa_cfw_f,$chain_r);

GetOptions(
		"i:s" => \$file,
		"o:s" => \$out_dir,
		"n:s" => \$name,
		"d" => \$d_flag,
		"s" => \$sature_f,
		"m" => \$mut_f,
		"mf" => \$mini_abund_mut,
		"c:s" => \$ref_s,
		"1:i" => \$one,
		"w:i" => \$window,
		"r:s" => \$ref_id,
		"nf"=> \$non_func_f,
		"lf"=> \$cdr3_len_f,
		"cf"=> \$aa_cfw_f,
		"t:s"=> \$chain_r
		);

die `pod2text $0` if (!$file && !$out_dir && !$name && !$ref_id);
$one = 3 unless(defined($one));
$window = 1000 unless(defined($window));
$mini_abund_mut = 10 unless(defined($mini_abund_mut));

my (%vusage, %jusage, %vj_pairing, %cdr3_length , %cdr3_freq , %clone_freq , %cdr3_aa, %clone_aa, %del_len, %ins_len,%vdj_len); # stat and output a file
my ($all_n, $cdr_aa_n , $mut_b_n , $all_b_n , $seq_mut_n) = (0,0,0 ,0 ,0,0);
my (%v_compo , %j_compo);
my $all_n_na = 0;



my %ref_sort;
# need to analyze the gene's base compostion 

if(defined($ref_s))
{
	open R, "$ref_s" or die;
	while(<R>)
	{
		chomp;
		s/^>//;
		my $id = $_;
		chomp(my $seq = <R>);
		my @s = split //, $seq;
		my $f=1;
		for(my $i=1 ; $i<=length($seq) ; $i++){
			if($s[$i-1] ne "-"){
				$ref_sort{$id}{$f} = $i;
				$f++;
			}
		}
	}
	close R;
}
#print Dumper(\%ref_sort);



#----------------------		read the file	------------------------
if($file=~/\.gz/){open I, "gzip -dc $file|" or die;}
else{open I, "$file" or die;}
$out_dir =~ s/\/$//;
my %flag_all;

<I>;
while(<I>)
{
	chomp;
	my @F = split;
	
	#	my $chain_flag = 0;
	if(defined($chain_r))# multiple chains analyzed together
	{
		my $chain_flag = 0;
		next if($F[2] eq "NA" || $F[4] eq "NA");# some CDR3 found by conserve rgion, without V or J
		if($chain_r eq "IGKL"){
			$chain_flag = 1 if(($F[2] =~ /^IGK/ || $F[2] =~ /^IGL/) && ($F[4] =~ /^IGK/ || $F[4] =~ /^IGL/));
		}else{
			$chain_flag = 1 if($F[2] =~ /^$chain_r/ && $F[4] =~ /^$chain_r/);
		}
		next if($chain_flag == 0);
	}
	#	next if($chain_flag == 0);

	$flag_all{$F[1]}++;
#	next if($F[1] eq "non-function");# sequence with peseudogene was filtered
#	next if($F[2] =~ /_unF/ or $F[4] =~ /_unF/);

	# filtration
	if($non_func_f){
		next if($F[1] =~ /out-of-frame/);
	}
	if($cdr3_len_f){
		next if(length($F[8])<3);
	}
	if($aa_cfw_f){
		my $f = 0;
		$f++ if($F[8] =~ /^C/);
		$f++ if($F[8] =~ /F$/ || $F[8] =~ /W$/);
		next if($f<2);
	}



        # Clone/CDR3 stat
        $cdr3_length{$F[6]-$F[5]+1}++;
        $cdr3_freq{$F[7]}++;
        if($F[8] ne "NA"){
                $cdr3_aa{$F[8]}++;
                $cdr_aa_n++;
        }
        $clone_freq{$F[17]}++;
        $clone_aa{$F[18]}++ if($F[18] ne "NA");
	
	$all_n_na++;
	next if($F[2] eq "NA" || $F[4] eq "NA");# some CDR3 found by conserve rgion, without V or J

	$all_n++;

	# V/J v-J usage
	my $vgene = (split /\*/,$F[2])[0];
	my $jgene = (split /\*/,$F[4])[0];
	$vusage{$vgene}++;
	$jusage{$jgene}++;
	$vj_pairing{"$vgene\t$jgene"}++;


	# deletion/insertion length stat
	if($F[9] eq "0"){$del_len{"V3"}{0}++;}
	else{$del_len{"V3"}{length($F[9])}++;}
	if($F[12] eq "0"){$del_len{"J5"}{0}++;}
	else{$del_len{"J5"}{length($F[12])}++;}

	if($d_flag && $F[29] ne "NA"){
		if($F[10] eq "0"){$del_len{"D5"}{0}++;}
		else{$del_len{"D5"}{length($F[10])}++;}
		if($F[11] eq "0"){$del_len{"D3"}{0}++;}
		else{$del_len{"D3"}{length($F[11])}++;}
		if($F[13] eq "0"){$ins_len{"VD"}{0}++}
		else{$ins_len{"VD"}{length($F[13])}++;}
		if($F[14] eq "0"){$ins_len{"DJ"}{0}++;}
		else{$ins_len{"DJ"}{length($F[14])}++;}
	}elsif(!$d_flag){
		if($F[15] eq "0"){$ins_len{"VJ"}{0}++;}
		else{$ins_len{"VJ"}{length($F[15])}++;}
	}

	# V,D,J gene length of CDR3 region
	my $v_len_cdr3 = $F[24]-$F[5]+1;
	my $j_len_cdr3 = $F[6]-$F[43]+1;
	$vdj_len{"V"}{$v_len_cdr3}++ if($v_len_cdr3 >= 0);
	$vdj_len{"J"}{$j_len_cdr3}++ if($j_len_cdr3 >= 0);
	if($d_flag && $F[29] ne "NA"){
		if($F[34]>$F[33]){$vdj_len{"D"}{$F[34]-$F[33]+1}++;}
		else{$vdj_len{"D"}{$F[33]-$F[34]+1}++;}
	}

	# calculate the hyper-mutation
#	if($mut_f)
#	{
#		my $v_mut = $F[22]-$F[23]+1;
#		my $v_b_n = $F[21]-$F[23]+1;
#		my $j_mut = $F[42]-(length($F[17])-$F[44]);
##		my $j_b_n = $F[41]-(length($F[17])-$F[44]);
#		$mut_b_n = $mut_b_n + $v_mut + $j_mut;
#		$all_b_n = $all_b_n + $v_b_n + $j_b_n;
#		$seq_mut_n++ if($v_mut + $j_mut>0);
#	}

	# analyze the gene's nucleotide composition
	if(defined($ref_s))
	{
		my $i;
		if($F[2] ne "NA")
		{
		my $v_s = substr($F[17],$F[23]-1,$F[24]-$F[23]+1);
		my @v_s_a = split //, $v_s;
		for($i=$F[25]; $i<=$F[26]; $i++){
			$v_compo{$ref_sort{$F[19]}{$i}}{$v_s_a[$i-$F[25]]}++;
		}
		}
		if($F[4] ne "NA")
		{
		my $j_s = substr($F[17],$F[43]-1,$F[44]-$F[43]+1);
		my @j_s_a = split //, $j_s;
		for($i=$F[45]; $i<=$F[46]; $i++){
			$j_compo{$ref_sort{$F[39]}{$i}}{$j_s_a[$i-$F[45]]}++;
		}
		}
	}
}
close I;



my $s_cdr3_aa = &output_sort("$out_dir/${name}_CDR3_AA.frequency.gz", 0, 1, %cdr3_aa);
#---------------        CDR3 aa frequency statistics    ------------------#
open I, "gzip -dc $out_dir/${name}_CDR3_AA.frequency.gz|" or die;
open S, ">$out_dir/${name}_CDR3_AA_section.stat" or die;
#open O, ">$out_dir/${name}_CDR3_AA_window.stat" or die;

my %uniq_sort;
#my %window_sort;
my $uniq_f = 0;

my ($s1, $s2 ,$s3 , $s4) = ("top100","top100-1000","top1000-1E4","more1E4");
for($s1, $s2 ,$s3 , $s4)
{
	$uniq_sort{$_} = [(0,0)];
}

my $cdr_aa_n_new = 0;
my %cdr3_aa_new;# used for saturation analysis
my %cdr3_aa_for_mut;

for (sort {$cdr3_aa{$b}<=>$cdr3_aa{$a}} keys %cdr3_aa)
{
	$uniq_f++;
	if($cdr3_aa{$_}>=$one)# less than $one, may many sequencing error
	{
		$cdr3_aa_new{$_} = $cdr3_aa{$_};
		$cdr_aa_n_new += $cdr3_aa{$_};
	}
	if($mut_f && $cdr3_aa{$_}>=$mini_abund_mut)# less than $mini_abund_mut, may many sequencing error
	{
		$cdr3_aa_for_mut{$_} = $cdr3_aa{$_};
	}

	my $rate = $cdr3_aa{$_}/$cdr_aa_n*100;
	if($uniq_f <= 100){
		$uniq_sort{$s1}->[0] = $rate;
		$uniq_sort{$s1}->[1] += $rate;
	}elsif($uniq_f <= 1000){
		$uniq_sort{$s2}->[0] = $rate;
		$uniq_sort{$s2}->[1] += $rate;
	}elsif($uniq_f <= 1E4){
		$uniq_sort{$s3}->[0] = $rate;
		$uniq_sort{$s3}->[1] += $rate;
	}else{
		$uniq_sort{$s4}->[0] = $rate;
		$uniq_sort{$s4}->[1] += $rate;
	}

#	my $section = int($cdr3_aa{$_}/$window);
#	$window_sort{$section}->[0]++;
#	$window_sort{$section}->[1] += $rate;

}

$uniq_sort{$s1}->[0]=sprintf("%1.0e",$uniq_sort{$s1}->[0]);
$uniq_sort{$s2}->[0]=sprintf("%1.0e",$uniq_sort{$s2}->[0]);
$uniq_sort{$s3}->[0]=sprintf("%1.0e",$uniq_sort{$s3}->[0]);
$uniq_sort{$s4}->[0]=sprintf("%1.0e",$uniq_sort{$s4}->[0]);
print S "top100\t>=$uniq_sort{$s1}->[0]%\t$uniq_sort{$s1}->[1]\n";
print S "100-1E3\t$uniq_sort{$s2}->[0]%-$uniq_sort{$s1}->[0]%\t$uniq_sort{$s2}->[1]\n";
print S "1E3-1E4\t$uniq_sort{$s3}->[0]%-$uniq_sort{$s2}->[0]%\t$uniq_sort{$s3}->[1]\n";
print S "more1E4\t>=$uniq_sort{$s4}->[0]%\t$uniq_sort{$s4}->[1]\n";


#-----=-------    Saturation curve(rarefraction curve)  ---------------------#
if($sature_f)
{
        open O, ">$out_dir/${name}_rarefraction_curve.txt" or die;
        print O "#Seq_num\tObserved\tChao1\tChao1_correct\n";
        my $firs = substr($cdr_aa_n_new,0,1);
        my $interval = $firs.(0 x (length($cdr_aa_n_new)-1));
        $interval = $interval/10;

        #1. 10 points
        for(my $i=1 ; $i<=10 ; $i++)
        {
                my ($S, $F1 , $F2) =  &get_random($cdr_aa_n_new , $interval*$i , $one , %cdr3_aa_new);#1. get the random sequences
                
		my ($chao1, $chao1_c) = ($S,$S);
                $chao1= $S + ($F1*$F1)/(2*$F2) if($F2);# Chao1 algorithm
                $chao1_c = $S + $F1*($F1-1)/(2*($F2+1)); # Chao1 corrected algorithm
                
		print O $interval*$i,"\t$S\t$chao1\t$chao1_c\n";
        }

        close O;
}

%cdr3_aa_new = ();

#-----------	Hyper-mutation analysis -----------------------------------#
my ($v_mut_all,$v_base_all,$d_mut_all,$d_base_all,$j_mut_all,$j_base_all) = (0,0,0,0,0,0);
my ($v_seq_mut,$d_seq_mut,$d_seq_all,$j_seq_mut,$vj_seq_all) = (0,0,0,0,0);
if($mut_f)
{
	if($file=~/\.gz/){open I, "gzip -dc $file|" or die;}
	else{open I, "$file" or die;}
	<I>;
	while(<I>)
	{
		chomp;
		my @F = split;
		next if($F[1] eq "non-function");# sequence with peseudogene was filtered
		next if($F[2] =~ /_unF/ or $F[4] =~ /_unF/);
		next if(!exists $cdr3_aa_for_mut{$F[8]});
		next if($F[2] eq "NA" || $F[4] eq "NA");

		$vj_seq_all++;
                my $v_mut = $F[22]-$F[23]+1;
                my $v_b_n = $F[21]-$F[23]+1;
		$v_seq_mut++ if($v_mut);
                my $j_mut = $F[42]-(length($F[17])-$F[44]);
                my $j_b_n = $F[41]-(length($F[17])-$F[44]);
		$j_seq_mut++ if($j_mut);
		
		my $d_mut = 0;
		my $d_b_n = 0;
		if($d_flag)# D gene
		{
			$d_mut = $F[32] if($F[32] ne "NA");
			$d_b_n = $F[31] if($F[31] ne "NA");
			$d_mut_all += $d_mut;
			$d_base_all += $d_b_n;
			$d_seq_mut++ if($d_mut);
			$d_seq_all++ if($F[31] ne "NA");
		}

		$v_mut_all += $v_mut;
		$v_base_all += $v_b_n;

		$j_mut_all += $j_mut;
		$j_base_all += $j_b_n;

#		print "$mut_b_n + $v_mut + $j_mut + $d_mut\n$all_b_n + $v_b_n + $j_b_n + $d_b_n\n";
                $mut_b_n = $mut_b_n + $v_mut + $j_mut + $d_mut;
                $all_b_n = $all_b_n + $v_b_n + $j_b_n + $d_b_n;
                $seq_mut_n++ if($v_mut + $j_mut + $d_mut >0);
#		print "$v_mut\t$v_b_n\t$j_mut\t$j_b_n\t$d_mut\t$d_b_n\n";
	}
	close I;
}




#for(sort {$a<=>$b} keys %window_sort)
#{
#	my $abund = $window*($_+1);
#	print O "$abund\t$window_sort{$_}->[0]\t$window_sort{$_}->[1]\n";
#}



#print Dumper(\%vdj_len);
#---------------	out put	-------------------------#
my $s_vusage = &output("$out_dir/${name}_V_gene.usage",$all_n, 1, %vusage);
my $s_jusage = &output("$out_dir/${name}_J_gene.usage",$all_n, 1, %jusage);
my $s_vjusage = &output("$out_dir/${name}_VJ_pairing.usage",$all_n ,1, %vj_pairing);

&output_sort_mini("$out_dir/${name}_CDR3_NT.length", $all_n_na, 0, %cdr3_length);
my $s_cdr3 = &output_sort("$out_dir/${name}_CDR3_NT.frequency.gz", $all_n_na, 1, %cdr3_freq);
my $s_clone = &output_sort("$out_dir/${name}_Clonotype_NT.frequency.gz",$all_n_na, 1, %clone_freq);
#my $s_cdr3_aa = &output_sort("$out_dir/${name}_CDR3_AA.frequency", 0, 1, %cdr3_aa);
my $s_clone_aa = &output_sort("$out_dir/${name}_Clonotype_AA.frequency.gz", 0, 1, %clone_aa);



#----- output the V,J gene's nucleotide composition
if(defined($ref_s))
{
	open V, ">$out_dir/${name}_V_NT_composition.txt" or die;
	open J, ">$out_dir/${name}_J_NT_composition.txt" or die;
	for my $p(sort {$a<=>$b} keys %v_compo){
		for (sort keys %{$v_compo{$p}}){
			print V "$p\t$_\t$v_compo{$p}{$_}\n";
		}
	}

	for my $p(sort {$a<=>$b} keys %j_compo){
		for (sort keys %{$j_compo{$p}}){
			print J "$p\t$_\t$j_compo{$p}{$_}\n";
		}
	}
}

#-------output the deletion length distribution
&output_len("$out_dir/${name}_VDJ_deletion_nt_len.dis", %del_len);
&output_len("$out_dir/${name}_VDJ_insertion_nt_len.dis", %ins_len);
&output_len("$out_dir/${name}_VDJ_length_inCDR3.dis", %vdj_len);



#--- print out the basical stat
my ($b_m_r,$s_m_r,$v_m_r,$v_s_r,$d_m_r,$d_s_r,$j_m_r,$j_s_r) = (0,0,0,0,0,0,0,0);
if($mut_f){
	open O, ">$out_dir/${name}_hypermutation.stat" or die;
	$v_m_r = $v_mut_all/$v_base_all*100 if($v_base_all);
	$v_s_r = $v_seq_mut/$vj_seq_all*100 if($vj_seq_all);

	$d_m_r = $d_mut_all/$d_base_all*100 if($d_base_all);
	$d_s_r = $d_seq_mut/$d_seq_all*100 if($d_seq_all);

	$j_m_r = $j_mut_all/$j_base_all*100 if($j_base_all);
	$j_s_r = $j_seq_mut/$vj_seq_all*100 if($vj_seq_all);

	$b_m_r = $mut_b_n/$all_b_n*100 if($all_b_n);
	$s_m_r = $seq_mut_n/$vj_seq_all*100 if($vj_seq_all);
	
	print O "V_gene $v_m_r\t$v_s_r\n";
	print O "D_gene $d_m_r\t$d_s_r\n";
	print O "J_gene $j_m_r\t$j_s_r\n";
	print O "Overall $b_m_r\t$s_m_r\n";
	close O;
}

my $tmp_all = 0;
for("in-frame","out-of-frame(stop_codon)","out-of-frame(CDR3_length)","non-function")
{
	$flag_all{$_} = 0 if(!exists $flag_all{$_});
	$tmp_all += $flag_all{$_};
}
my $r = 0;
$r = sprintf("%0.2f",$flag_all{'in-frame'}/$tmp_all*100) if($tmp_all);
print "in-frame: $flag_all{'in-frame'}\t$r\n";
$r=  sprintf("%0.2f",$flag_all{'out-of-frame(stop_codon)'}/$tmp_all*100) if($tmp_all);
print "out-of-frame(stop_codon): $flag_all{'out-of-frame(stop_codon)'}\t$r\n";
$r=  sprintf("%0.2f",$flag_all{'out-of-frame(CDR3_length)'}/$tmp_all*100) if($tmp_all);
print "out-of-frame(CDR3_length): $flag_all{'out-of-frame(CDR3_length)'}\t$r\n";
$r=  sprintf("%0.2f",$flag_all{'non-function'}/$tmp_all*100) if($tmp_all);
print "non-function: $flag_all{'non-function'}\t$r\n";


open I, "$ref_id" or die;
my %ref_id_gene;
while(<I>)
{
	chomp;
	my @a = split /:/,$_;
	for(@a)
	{
		next if(/_unF/);
		my $g = (split /\*/,$_)[0];
		if(/V/){
			$ref_id_gene{"V"}{$g} = 1;
		}elsif(/J/){
			$ref_id_gene{"J"}{$g} = 1;
		}
	}
}
close I;
my ($ref_v_n , $ref_j_n) = (scalar keys %{$ref_id_gene{'V'}},scalar keys %{$ref_id_gene{'J'}});
my $ref_vj_n = $ref_v_n*$ref_j_n;
$ref_v_n = scalar(keys %vusage)/$ref_v_n*100;$ref_v_n=sprintf("%0.2f",$ref_v_n);
$ref_j_n = scalar(keys %jusage)/$ref_j_n*100;$ref_j_n=sprintf("%0.2f",$ref_j_n);
$ref_vj_n = scalar(keys %vj_pairing)/$ref_vj_n*100;$ref_vj_n=sprintf("%0.2f",$ref_vj_n);
print "V_gene_used: ",scalar(keys %vusage),"\t$ref_v_n\n";
print "J_gene_used: ",scalar(keys %jusage),"\t$ref_j_n\n";
print "V-J_pairing: ",scalar(keys %vj_pairing),"\t$ref_vj_n\n";
print "Uniq_number(seq_nt,seq_aa): ",scalar(keys %clone_freq),"\t",scalar(keys %clone_aa),"\n";
print "Uniq_number(cdr3_nt,cdr3_aa): ",scalar(keys %cdr3_freq),"\t",scalar(keys %cdr3_aa),"\n";;
($s_clone,$s_clone_aa,$s_cdr3,$s_cdr3_aa) = (sprintf("%0.2f",$s_clone),sprintf("%0.2f",$s_clone_aa),sprintf("%0.2f",$s_cdr3),sprintf("%0.2f",$s_cdr3_aa));
print "Shannon_index(seq,seq_aa): $s_clone\t$s_clone_aa\n";
print "Shannon_index(cdr3_nt,cdr3_aa): $s_cdr3\t$s_cdr3_aa\n";
($s_vusage,$s_jusage,$s_vjusage) = (sprintf("%0.2f",$s_vusage),sprintf("%0.2f",$s_jusage),sprintf("%0.2f",$s_vjusage));
print "Shanono_index(V,J,V-J): $s_vusage\t$s_jusage\t$s_vjusage\n";
($b_m_r,$s_m_r) = (sprintf("%0.2f",$b_m_r),sprintf("%0.2f",$s_m_r));
print "Hyper-mutation(base_rate,seq_rate): $b_m_r\t$s_m_r\n";



#--------------------#
#  get random sequence
#--------------------#
sub get_random
{
	my ($sum , $num , $low , %h) = @_;
	my %random;
	my ($uniq_num ,$F1 , $F2) = (0,0,0);
	# get the sequence id selected at random
	if($num >= $sum) # if the $sum is less than $num, all $sum id is selected
	{
		$random{$_} = 1 for(1..$sum);
	}
	else{
		while(1)
		{
        		my $f = int rand($sum);
        		$random{($f+1)}=1;
        		last if(scalar keys %random == $num);
		}
	}

	# get the selected sequences' uniq number , one , two
	my $flag = 1;
	for(keys %h)
	{
        	my $f = 0;
        	for(my $i=$flag; $i<$flag+$h{$_} ; $i++){
                	$f++ if(exists $random{$i});
        	}
        	$flag += $h{$_};
		if($f){
			$uniq_num++;
		}
		if($f==$low){
			$F1++;
		}
		if($f==$low+1){
			$F2++;
		}
	}
	return($uniq_num,$F1,$F2);

}



#----------------#
# 
#----------------

sub output_len
{
	my ($out, %h) = @_;
	open O, ">$out" or die;
	for(my $j=0;$j<=30;$j++){print O "\t$j";}
	print O"\n";

	for my $g(sort {$b cmp $a} keys %h)
	{
		my $sum = 0;
		print O "$g";
		for (keys %{$h{$g}}){
			$sum += $h{$g}{$_};
		}
		for (my $i=0 ; $i<= 30 ; $i++){
			$h{$g}{$i} = 0 unless(exists($h{$g}{$i}));
			my $rate = $h{$g}{$i}/$sum*100;
			print O "\t$rate";
		}
		print O "\n";
	}
	close O;
}


#------------#
# output sub
#------------#
# sort by the value
# calcultate the Shannon Index
sub output_sort{
	my ($outfile , $sum, $shannon_f, %hash) = @_;
	my $shannon = 0;
	open O, "|gzip>$outfile" or die;
	if($sum==0){
		$sum += $hash{$_} for(keys %hash);
	}
	for (sort {$hash{$b} <=> $hash{$a}} keys %hash){
		my $rate = $hash{$_}/$sum;
		if($shannon_f){
			$shannon = $shannon-$rate*(log($rate)/log(2));
		}
		print O "$_\t$hash{$_}\t",$rate*100,"\n";
	}
	close O;
	return $shannon if($shannon_f);
}

#sort by the key
sub output{
	my ($outfile , $sum,$shannon_f, %hash) = @_;
	my $shannon = 0;
	open O, ">$outfile" or die;
	for (sort keys %hash){
		my $rate = $hash{$_}/$sum;
		if($shannon_f){
			$shannon = $shannon-$rate*(log($rate)/log(2));
		}
		print O "$_\t$hash{$_}\t",$rate*100,"\n";
	}
	close O;
	return $shannon if($shannon_f);
}

# sort by the key
sub output_sort_mini
{
	my ($outfile , $sum,$shannon_f, %hash) = @_;
	my $shannon = 0;
	open O, ">$outfile" or die;
	for (sort {$a<=>$b} keys %hash){
		my $rate = $hash{$_}/$sum;
		if($shannon_f){
			$shannon = $shannon-$rate*(log($rate)/log(2));
		}
		print O "$_\t$hash{$_}\t",$rate*100,"\n";
	}
	close O;
	return $shannon if($shannon_f);

}
