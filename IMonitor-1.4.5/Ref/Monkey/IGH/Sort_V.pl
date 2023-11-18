#!/usr/bin/perl -w
use strict;

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


open I, "IGHV.MacVJ.raw.fa" or die;
my %All;
my $max = 0;
while(<I>)
{
	chomp;
	my $id = $_;
	#print "$_\n";
	chomp(my $seq = <I>);
	my $len = length($seq);
	my $flag = 0;
	for(my $i=$len;$i>=9;$i--){
		my $A1 = substr($seq,$i-9,3);
		my $A3 = substr($seq,$i-3,3);
		if($codon2aa{$A1} eq "Y" and $codon2aa{$A3} eq "C"){
			$All{"$id:$seq"} = $i;
			$max = $i if($i > $max);
			$flag = 1;
			last;
		}
	}
	if($flag == 0){
		print "Error: $id\n";
	}
}
close I;


for(keys %All)
{
	my ($id, $seq) = split /:/,$_;
	if($All{$_} == $max){
		print "$id\n$seq\n";
	}else{
		my $add = $max-$All{$_};
		my $dot = "." x $add;
		my $seq = $dot.$seq;
		print "$id\n$seq\n";
	}
}
print "$max\n";
