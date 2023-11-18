#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;


=head1 Usage:

	perl Use_primer_find_pcr_seq.pl <*.primer> <*.fa> <[V/J]> <out>

=head1 Introduction:
	
	Align the primers to the *.fa file and find the PCR templates with some mismatch, then print out the sequence "primer + amplified sequence". One template may find multiple corresponding primers. 
	for example, primer:RCT, tmplate seq: AAATACTGGGC. 
	then print out 2 sequences: ACTGGGC GCTGGGC

=head1 Input file format:
	
	<*.primer>: pirmer_id       primer_seq e.g. IGHV1-46/2-F2   GTCACCATGACCAGGGACACG
	<*.fa>: normal fasta format from IMGT, must be sorted e.g. :
		>J00256|IGHJ3*01|Homo sapiens|F|J-REGION|1537..1586|50 nt|2| | | | |50+0=50| | |
		-----TGATGC---TT------TTGATGTCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG
	NOte: the whole primer must be included in the template sequence. or it will not be considered. 
	     e.g.: template seq: AAAACCgtac, primter seq: GTACG, in this condition, no result given


=head1 print out file format:

	>gene_type_name:numbering:primer_start_position:out_seq_length:primer_len:primer_id:mismatch:t/p(t:template,p:primer)
	out_sequence
	e.g.
	>IGHV5-a*02:1:237:83:19:IGH3:3:p
	AGCTGACAAGTCCATCAGCACTGCCTACCTGCAGTGGAGCAGCCTGAAGGC-TCGGACACCGCCATGTATTACTGTGCGAGACA
	
=cut

unless (@ARGV == 4)
{
	die `pod2text $0`;
}
my ($in_file1 , $in_file2 , $gene , $out) = @ARGV;
open IN1 , "$in_file1" or die;
open IN2 , "$in_file2" or die;
open OUT , ">$out" or die;

my $mismatch_num = 2;
my $V_restrict = 100;
$V_restrict = 0 if($gene eq "J");

my %degeneration = ("A"=>"A","C"=>"C","T"=>"T","G"=>"G",
		    "M"=>"AC","R"=>"AG","W"=>"AT","S"=>"CG","Y"=>"CT","K"=>"GT",
		    "V"=>"ACG","H"=>"ACT","D"=>"AGT","B"=>"CGT","N"=>"ACGT");

my %primer;

#----------------------------- read Primer file
# $in_file1 format:
#	pirmer_id	primer_seq. e.g. IGHV1-46/2-F2   GTCACCATGACCAGGGACACG

while(<IN1>)
{
	chomp;
	my ($id,$seq) = split;
	next if($gene eq "V" and $id =~ /J/);
	next if($gene eq "J" and $id =~ /V/);
	$seq =~ tr/acgt/ACGT/;

#	if($gene eq "J")# for J gene , change primer sequence into complementary strand
#	{
#		$seq =~ tr/ACGT/TGCA/;
#		$seq = reverse $seq;
#	}

	if($seq =~ /[^ACGT]/)# has degenerate base
	{
		my %mulit;
		my @temp = split // , $seq;
		for my $base(@temp)
		{
			if(scalar keys %mulit == 0)
			{
				for (split //,$degeneration{$base})
				{
					$mulit{$_} = 1;
				}
			}
			else
			{
				for my $part_seq (keys %mulit)
				{
					$mulit{"$part_seq$_"} = 1 for (split //,$degeneration{$base});
					delete $mulit{$part_seq};
				}
			}
		}

		for (keys %mulit)
		{
			$primer{length($seq)}{$_} = $id ;
		}
	}
	else
	{
		$primer{length($seq)}{$seq} = $id;
	}
}
close IN1;

#print Dumper(\%primer);

#-------------	mismatch tolerant	---------------------
#for (my $i=0 ; $i<$mismatch_num ; $i++)
#{
#	&one_mismatch_tolerant(\%primer);
#}
#print Dumper(\%primer);


my @sort_seq_out;# store all sequences needed to print out and then sort them 
#my $max_len = 0;
my ($min_start , $max_loc ) = (1000 , 0);

$/=">";
<IN2>;
while(<IN2>)
{
	chomp;
	my @line = split /\n+/,$_;
#	my $id = shift @line;
	my $id = (split /\|/,$line[0])[1];
	if((split /\|/, $line[0])[3] ne "F" and (split /\|/, $line[0])[3] ne "(F)" and (split /\|/, $line[0])[3] ne "[F]"){$id.="_unF";}# # for pseudogene , add a tag

	shift @line;
	my $raw_seq="@line";
	$raw_seq=~s/\s+//g;

#	next if(/^$/);
#	my $id = (split /\|/,$_)[1];
#	chomp(my $raw_seq = <IN2>);

	$raw_seq =~ tr/acgt/ACGT/;
	if($gene eq "J") # for J gene , change the sequence into its complementary strand
	{
		$raw_seq =~ tr/ACGT/TGCA/;
		$raw_seq = reverse $raw_seq;
	}

	$max_loc = length($raw_seq) if($max_loc < length($raw_seq));# for next sorting

	$raw_seq =~s/-/./g;
	my $rm_dot_seq = $raw_seq; $rm_dot_seq =~ s/\.//g;# get the sequence removed dots

	my %coordinates; # get the coordinates between raw_seq and removed dots' sequence
	$V_restrict = &P_to_P_cor($raw_seq , $rm_dot_seq , \%coordinates);
#	print Dumper(\%coordinates);
	
	# find the primer
	my %new_seq_out;

	for my $len(keys %primer)
	{
		for (my $i=$V_restrict ; $i<=length($rm_dot_seq)-$len ; $i++)# start to find template after $V_restrict postion
		{
			my $subseq = substr($rm_dot_seq , $i , $len);
			for (sort keys %{$primer{$len}})
			{
				my $flag = &two_seq_align($subseq , $_);# compare two subsequences
				if($flag <= $mismatch_num)
				{
					my $new_seq_len = length($rm_dot_seq) - $i;
					my $new_seq_1 = substr($raw_seq , $coordinates{$i});# template amplify

					# get the sequence: primer + amplify
					my $new_seq_2 = &primer_sorted($_ ,$new_seq_1) ; # primer sequence sorted based on corresponding region sequence of template


					my $primer_start_pos = $coordinates{$i}+1;
					$min_start = $primer_start_pos if($min_start > $primer_start_pos);# get the min start postion for next sorting

					
					$new_seq_out{scalar keys %new_seq_out}{"$id:".(scalar keys %new_seq_out).":".$primer_start_pos.":$new_seq_len:$len:$primer{$len}{$_}:0:t:$subseq:$subseq"} = $new_seq_1;
					if($flag==0){
						$new_seq_out{scalar keys %new_seq_out}{"$id:".(scalar keys %new_seq_out).":".$primer_start_pos.":$new_seq_len:$len:$primer{$len}{$_}:$flag:p:$subseq:$subseq"} = $new_seq_2;
					}
					else{
						$new_seq_out{scalar keys %new_seq_out}{"$id:".(scalar keys %new_seq_out).":".$primer_start_pos.":$new_seq_len:$len:$primer{$len}{$_}:$flag:p:$_:$subseq"} = $new_seq_2;
					}

				}
			}
		}
	}
	
	# store in one array
	for my $numbering(sort {$a <=> $b} keys %new_seq_out)
	{
		for (keys %{$new_seq_out{$numbering}})
		{
			push @sort_seq_out , "$_\t$new_seq_out{$numbering}{$_}";
		}
	}
}

$/="\n";



#  sort the sequences and print out;
for (my $i=0 ; $i<=$#sort_seq_out ; $i++)
{
	my ($seq_id , $new_seq)= split /\t/,$sort_seq_out[$i];
	my $p_cut_pos = (split /:/ , $seq_id)[2];
		
	$new_seq = ("-" x ($p_cut_pos - $min_start)).$new_seq;# sorting: 5' region 
	$new_seq .= "-" x ($max_loc - length($new_seq) - $min_start + 1);# sorting: 3' region

	if($gene eq "J")# for J, change the sequence into complementary strand
	{
		$new_seq = reverse $new_seq;
		$new_seq =~ tr/ACGT/TGCA/;
		$p_cut_pos = length($new_seq)-$p_cut_pos+1; # primer cut position changed because of strand chenged
		my @temp = split /:/ , $seq_id; $temp[2] = $p_cut_pos; $seq_id = join ":" , @temp;

	}
	$new_seq =~ s/\./-/g;
	print OUT ">$seq_id\n$new_seq\n";
}





#---------------------------------------------------------------------------------------#
#	the corresponding coordinates between sequence and sequences without dots	#
#---------------------------------------------------------------------------------------#

sub P_to_P_cor
{
	my ($seq1 , $seq2 , $cor) = @_;
	my $flag = 0;

	my @seq_1 = split // , $seq1;
	my @seq_2 = split // , $seq2;
	my $j=0;
	for (my $i=0 ; $i<=$#seq_2 ; $i++)
	{
		for(;$j<=$#seq_1 ;)
		{
			$flag = $i if($j eq $V_restrict);
			if($seq_1[$j] eq "." || $seq_1[$j] eq "-")
			{
				$j++;
				next;
			}
			$$cor{$i} = $j;
				
			$j++;
			last;
		}
	}
	return $flag;
}


#-----------------------------------------------------#
# compare two sequences and return mismatch number----#
#-----------------------------------------------------#
sub two_seq_align
{
	my ($seq1 , $seq2) = @_;
	my @seq_1 = split // , $seq1;
	my @seq_2 = split // , $seq2;
	my $flag = 0;
	for (my $i=0 ; $i<=$#seq_1 ; $i++)
	{
		$flag++ if($seq_1[$i] ne $seq_2[$i]);
		last if($flag > $mismatch_num);
	}
	return $flag;
}


#-------------------------------#
# primer sequence sorted	#
#-------------------------------#
sub primer_sorted
{
	my ($p , $s) = @_;
	my @s_all = split // , $s;
	my @p_all = split // , $p;
	for (@s_all)
	{
		if($_ ne "-" and $_ ne "."){
			$_ = $p_all[0];
			shift @p_all;
			last if(scalar @p_all == 0);
		}
	}

	$p = join "" , @s_all;
	return $p;
}
