#!/usr/bin/perl -w
use strict;

unless(@ARGV==2){
	print "perl $0 <in><out><stat>\n";
	print "Instruction for translated amino acid:\n";
	print "-,codon with 'N';*,stop codon;x,codon less than 3 nt at the end\n\n";
	exit;
}



open O, ">$ARGV[1]" or die;


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


my %conserve_j;
#$conserve_j{1} = "W";
$conserve_j{2} = "G";
$conserve_j{4} = "G";
my %conserve_v;
$conserve_v{1}= "Y";
#$conserve_v{2}= "Y";
$conserve_v{3}= "C";

open I, "gzip -dc $ARGV[0]|" or die;

my $find_cdr3=0;
my ($sum , $with_v , $with_j, $without_vj) = (0,0,0,0);
while(<I>)
{
	chomp;
	my $flag = "in-frame";
	my $strand = "Plus";
	my @F = split;
	my ($cdr3_start, $cdr3_end) = (0,0);
	
	$sum++;
	if($F[0] ne "NA" && $F[20] ne "NA"){
		next;
	}

	my $q_len = length($F[31]);
	my $tmp;
	if($F[0] ne "NA" && $F[20] eq "NA")# only V was mapped 
	{
		if(($F[6]-$F[7])>0)# change into plus strand
		{
			$strand="Minus";
	                $F[31]=reverse $F[31];
	                $F[31]=~tr/ATCGatcg/TAGCtagc/;
	                $F[4]=$q_len-$F[4]+1;
	                $F[5]=$q_len-$F[5]+1;

	                $tmp=$F[4];$F[4]=$F[5];$F[5]=$tmp;
	                $tmp=$F[6];$F[6]=$F[7];$F[7]=$tmp;
		}
		my $cdr3_ref_start = (split /\./,$F[0])[4];
		$cdr3_start = $F[5]-($F[7]-$cdr3_ref_start);
		$cdr3_end = &search_cdr3_by_aa($F[31],"J");
		$with_v++ if($cdr3_end && $cdr3_end-$cdr3_start+1 > 0);

	}elsif($F[0] eq "NA" && $F[20] ne "NA")# only J was mapped
	{
		if(($F[26]-$F[27])>0)# change into plus strand
		{
			$strand="Minus";
       		        $F[31]=reverse $F[31];
	                $F[31]=~tr/ATCGatcg/TAGCtagc/;
	                $F[24]=$q_len-$F[24]+1;
	                $F[25]=$q_len-$F[25]+1;

	                $tmp=$F[24];$F[24]=$F[25];$F[25]=$tmp;
			$tmp=$F[26];$F[26]=$F[27];$F[27]=$tmp;
		}
		my $cdr3_ref_end = (split /\./,$F[20])[4];
		$cdr3_end = $F[24] + ($cdr3_ref_end-$F[26]);
		$cdr3_start = &search_cdr3_by_aa($F[31],"V");
		$with_j++ if($cdr3_start && $cdr3_end-$cdr3_start+1 > 0);

	}elsif($F[0] ne "NA" && $F[20] ne "NA")# Both V and J was not mapped
	{
		$cdr3_start = &search_cdr3_by_aa($F[31],"V");
		if($cdr3_start == 0){
			$strand="Minus";
			$F[31]=reverse $F[31];
			$F[31]=~tr/ATCGatcg/TAGCtagc/;
			$cdr3_start = &search_cdr3_by_aa($F[31],"V");
			$cdr3_end = &search_cdr3_by_aa($F[31],"J");
		}
		$without_vj++ if($cdr3_start && $cdr3_end && $cdr3_end-$cdr3_start+1 > 0);
	}

	if($cdr3_start==0 || $cdr3_end==0){
		next;
	}

	#  CDR3 filter && translate into AA
	my $cdr3length=$cdr3_end-$cdr3_start+1;
	if($cdr3length<=0){
		next;
	}
	
	my $V_phase = ($cdr3_start-1)%3;
	my $nt_for_clonotypeAA = $F[31];
	my $aa_for_clonotypeAA = "NA";
	my $aa = "NA";
	my $ORFdna = substr($F[31],$cdr3_start-1,$cdr3length);
	if($cdr3length%3!=0){
		$flag = "out-of-frame(CDR3_length)";
	}else{
		$aa=&dna2aa($ORFdna,0);
		$aa_for_clonotypeAA=&dna2aa($nt_for_clonotypeAA,$V_phase);
		$flag = "out-of-frame(stop_codon)" if($aa_for_clonotypeAA=~/\*/);
	}

	# out put
	$find_cdr3++;
	my $id = $F[30];
	pop @F;
	pop @F;
	my $new_line = join "\t" , @F;

	print O	"$id\t$flag\t$F[0]\t$F[10]\t$F[20]\t$cdr3_start\t$cdr3_end\t$ORFdna\t$aa\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$strand\t$nt_for_clonotypeAA\t$aa_for_clonotypeAA\t$new_line\n";
}
close I;

print "All_seq: $sum\n";
print "Find_CDR3_byconserve: $find_cdr3\t",$find_cdr3/$sum*100,"\n";
print "Find_CDR3_withV(byconserve): $with_v\t",$with_v/$sum*100,"\n";
print "Find_CDR3_withJ(byconserve): $with_j\t",$with_j/$sum*100,"\n";
print "Find_CDR3_withoutVJ(byconserve): $without_vj\t",$without_vj/$sum*100,"\n";

#-------------#
# search CDR3 by conserve aa
#---------------------------#

sub search_cdr3_by_aa
{
	my ($subseq , $gene) = @_;
	if($gene eq "V")
	{
		for(my $i=0 ; $i<=length($subseq)-9; $i++)
		{
			my $s1 = substr($subseq, $i, 3);
			my $s2 = substr($subseq, $i+3, 3);
			my $s4 = substr($subseq, $i+6, 3);
			my $flag = 0;
			for(($s1,$s2,$s4)){
				$codon2aa{$_} = "-" if(!exists $codon2aa{$_});
			}
			$flag++ if($codon2aa{$s1} eq $conserve_v{1});
#			$flag++ if($codon2aa{$s2} eq $conserve_v{2});
			$flag++ if($codon2aa{$s4} eq $conserve_v{3});
			if($flag==3){
				return($i+7);
			}
		}
	}else{
                for(my $i=0 ; $i<=length($subseq)-12; $i++)
                {
                        my $s1 = substr($subseq, $i, 3);
                        my $s2 = substr($subseq, $i+3, 3);
                        my $s4 = substr($subseq, $i+9, 3);
                        my $flag = 0;
			for(($s1,$s2,$s4)){
				$codon2aa{$_} = "x" if(!exists $codon2aa{$_});
			}
                        $flag++ if($codon2aa{$s1} eq "F" || $codon2aa{$s1} eq "W");
                        $flag++ if($codon2aa{$s2} eq $conserve_j{2});
                        $flag++ if($codon2aa{$s4} eq $conserve_j{4});
                        if($flag==3){
                                return($i+12);
                        }
                }
	
	}
}


# -------------#
# translate
#--------------#

sub dna2aa{
        my $seq=shift;
#number of nucleotide before condon in the $seq
        my $phase=shift;

        $phase==1 and $seq=~s/^\w//;
        $phase==2 and $seq=~s/^\w\w//;
        my $aa="";
        for (my $start=0;$start<length($seq);$start+=3){
                my $codon=uc(substr($seq,$start,3));
                !exists $codon2aa{"$codon"} and $codon2aa{"$codon"}="x";
                $aa.=$codon2aa{"$codon"};
        }
        return $aa;

}

