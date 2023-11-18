#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;



=head1 Usage:

	perl Ref_seq_integration.pl <*V.fa> [<*D.fa>] <*J.fa> <CDR3-region> <type> <out_dir>

=head1 Input files:
	
	<*V.fa> <*D.fa> <*J.fa> : all sorted and the format as follow:
	e.g.
	<*V.fa>: 
	>IGHV1-18*01:1:223:98:22:IGHV1-24-F2:0:t
	------AGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGA
	<*D.fa>:
	>J00235|IGHD2-21*01|Homo
	-AGCATATTGTGGTG---------GTGAT---TGCTATTCC
	<*J.fa>:
	>IGHJ3*01:1:64:50:22:IGHJ-CONS-2:0:t
	-----TGATGC---TT------TTGATGTCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG
	<CDR3-region>:

	<type>: e.g. IGH,IGKL,TRB,TRA..

=head1 Output:

	create 7(without D gene) or 9 files in the directory <out_dir>.the files' names as follow:
	*.[VDJ].converted.fa
	*.[VDJ].alignment
	*.converted.fa
	*.converted.SingleLine
	*.backup

	out new id : sub_family.same_seq_num.real_seq_len.primer_len.CDR3_real_start_pos.numbering
		     TRBV6-6*04.1.103.23.89.186
	Note:	for the part sequence without CDR3, the give the align potential position. 
		e.g.  ACGTGCG-----, if the aligned CDR3 position was 10, the give the potential position 10
=cut

die `pod2text $0` if(@ARGV != 5 and @ARGV != 6);

my ($V_file , $D_file , $J_file , $cdr3_region , $type , $out_dir);
($V_file , $D_file , $J_file , $cdr3_region , $type , $out_dir) = @ARGV if(@ARGV == 6);# has D gene
($V_file , $J_file , $cdr3_region , $type , $out_dir) = @ARGV if(@ARGV == 5);# withtou D gene
$out_dir =~ s/\/$//;

my $flag = 0;
$flag = 1 if(@ARGV == 6); # the flag with D gene

my ($CDR3_V , $CDR3_J) = (split /J/ , $cdr3_region); # get CDR3 region on template
$CDR3_V =~ s/^V//;

# if out directory exists *.backup, *.converted.SingleLine , *.converted.fa files, then remove them
unlink "$out_dir/$type.backup";
unlink "$out_dir/$type.converted.SingleLine";
unlink "$out_dir/$type.converted.fa";

# print out files : *.backup, *.converted.SingleLine , *.converted.fa 
open O_B, ">>$out_dir/$type.backup" or die;
open O_S, ">>$out_dir/$type.converted.SingleLine" or die;
open O_C, ">>$out_dir/$type.converted.fa" or die;


my %V_ref_uniq; # key: combined id, value: sequence
my %D_ref_uniq; # key: combined id, value: sequence
my %J_ref_uniq; # key: combined id, value: sequence

&integration_sequences_VJ(\%V_ref_uniq , $V_file , $CDR3_V); # for V ,get the unique sequences and composite the id; print out *.backup files
#print Dumper(\%V_ref_uniq);
&Print_out(\%V_ref_uniq , "V");# print out *.alignment, *.converted.fa, *.converted.SingleLine
&integration_sequences_VJ(\%J_ref_uniq , $J_file , $CDR3_J); # for J ,get the unique sequences and composite the id; print out *.backup files
&Print_out(\%J_ref_uniq , "J");# print out *.alignment, *.converted.fa, *.converted.SingleLine
&integration_sequences_D(\%D_ref_uniq , $D_file) if($flag); # for J ,get the unique sequences and composite the id; print out *.backup files
&Print_out(\%D_ref_uniq , "D") if($flag);# print out *.alignment, *.converted.fa, *.converted.SingleLine

close O_C;
close O_B;
close O_S;

#-------------------------------------------------------#
# for VJ, get the unique sequences and composite the id	#
#-------------------------------------------------------#
sub integration_sequences_VJ
{
	my ($ref_uniq , $file , $CDR3_edge) = @_;
	open IN, "$file" or die;	
	my %all_seq;

	# read file and store in a hash
	my $numbering = 0;
	my %id_uniq;
	while(<IN>)
	{
		next if(/^$/);
		chomp;
		my $raw_id = $_;
		$raw_id =~ s/^>//;
		my ($id , $len , $primer_len , $prim_id , $mis , $tag ,$p_seq , $ref_seq) = (split /:/ , $raw_id)[0,3,4,5,6,7,8,9];
		($p_seq , $ref_seq) = ("N","N") if(!defined($p_seq));
		chomp(my $seq = <IN>);
		$seq =~ tr/acgt/ACGT/;
		# get CDR3 position 
		my $sub_seq = substr($seq , 0 , $CDR3_edge);
		my $non_base_num = $sub_seq =~ tr/-/-/;

		if($sub_seq =~ /-$/)# for the sequence with part seq.
		{
			$sub_seq =~ /(-+$)/;
			$non_base_num -= length($1);
		}
		my $CDR3_real_edge = $CDR3_edge - $non_base_num;
		
		# store unique sequence and composite sequences' id
		if(!exists $all_seq{"$primer_len:$p_seq:$ref_seq:$seq"})
		{
			$all_seq{"$primer_len:$p_seq:$ref_seq:$seq"} = [($id , $numbering , $len , $primer_len , $CDR3_real_edge)];
			$numbering++;
			$id_uniq{"$primer_len:$p_seq:$ref_seq:$seq"}{$id} = 1;
		}
		else
		{
			if(!exists $id_uniq{"$primer_len:$p_seq:$ref_seq:$seq"}{$id}) # combine the sequences' id
			{
				$all_seq{"$primer_len:$p_seq:$ref_seq:$seq"}->[0] .= ":$id";
				$id_uniq{"$primer_len:$p_seq:$ref_seq:$seq"}{$id} = 1;
			}
			$all_seq{"$primer_len:$p_seq:$ref_seq:$seq"}->[3] = $primer_len if($all_seq{"$primer_len:$p_seq:$ref_seq:$seq"}->[3]<$primer_len); # the longer primer was selected
		}

	}
	close IN;
	
	# conduct multiple type ids
	for my $seq1 ( sort keys %all_seq)
	{

		my ($p_l ,$p_s , $r_s , $s) = split /:/, $seq1;
		my @ids = split /:/ , $all_seq{$seq1}->[0];
		my $id_num = scalar @ids;# the number of types with same sequences
		my $sub_family; # gene's sub-family
		my %temp_id;
		$temp_id{(split /\*/ , $_)[0]} = 1 for (@ids);
	
		if($id_num == 1) # has only one gene type
		{
			$sub_family = $ids[0];
		}
		else # has multiple types
		{
			$sub_family = join "" , keys %temp_id;
		}

		# store combined information in a hash
		my $new_id = "$sub_family.$id_num.$all_seq{$seq1}->[2].$all_seq{$seq1}->[3].$all_seq{$seq1}->[4].$all_seq{$seq1}->[1]";
		$$ref_uniq{$new_id} = $s;#  store combined information in a hash

		# print out *.backup file
		print O_B "$all_seq{$seq1}->[0]\t$new_id\n";
		print "$new_id\t$p_l\t$p_s\t$r_s\n";
	}


}


#---------------------------------------------------------------#
#  for D , get the unique sequences and composite the id	#
#---------------------------------------------------------------#

sub integration_sequences_D
{

	my ($ref_uniq , $file) = @_;
	open IN, "$file" or die;	
	my %all_seq;

	# read file and store in a hash
	my $numbering = 0;
	
	while(<IN>)
	{
		next if(/^$/);
		chomp;
		my $raw_id = $_;
		$raw_id =~ s/^>//;
		my $id = (split /\|/ , $raw_id)[1];

		chomp(my $seq = <IN>);
		$seq =~ tr/acgt/ACGT/;
		# get sequence length
		my $len = $seq;
		$len =~ s/-//g;
		$len = length($len);

		# store unique sequence and composite sequences' id
		if(!exists $all_seq{$seq})
		{
			$all_seq{$seq} = [($id , $numbering , $len , 0 , 0)];
			$numbering++;
		}
		else
		{
			$all_seq{$seq}->[0] .= ":$id";
		}

	}
	close IN;
	
	# conduct multiple type ids
	for my $seq1 ( sort keys %all_seq)
	{
		my @ids = split /:/ , $all_seq{$seq1}->[0];
		my $id_num = scalar @ids;# the number of types with same sequences
		my $sub_family; # gene's sub-family
		my %temp_id;
		$temp_id{(split /\*/ , $_)[0]} = 1 for (@ids);
	
		if($id_num == 1) # has only one gene type
		{
			$sub_family = $ids[0];
		}
		else # has multiple types
		{
			$sub_family = join "" , keys %temp_id;
		}

		# store combined information in a hash
		my $new_id = "$sub_family.$id_num.$all_seq{$seq1}->[2].$all_seq{$seq1}->[3].$all_seq{$seq1}->[4].$all_seq{$seq1}->[1]";
		$$ref_uniq{$new_id} = $seq1;#  store combined information in a hash

		# print out *.backup file
		print O_B "$all_seq{$seq1}->[0]\t$new_id\n";
	}
}


#-----------------------------#
#	print out files	      #
#-----------------------------#
# print out :
# *.[VDJ].converted.fa
# *.[VDJ].alignment
# *.converted.fa
# *.converted.SingleLine

sub Print_out
{
	my $ref_uniq = shift;
	my $gene = shift;

	open C , ">$out_dir/$type.$gene.converted.fa" or die;# create file handle for *.converted.fa
	open A , ">$out_dir/$type.$gene.alignment" or die; # create file handle for *.alingment
	
	# print out .....
	for my $new_id (keys %$ref_uniq)
	{
		print A ">$new_id\n$$ref_uniq{$new_id}\n";# print out *.alingment

		my $new_seq = $$ref_uniq{$new_id};
		$new_seq =~ s/-//g;# remove bars
		print C ">$new_id\n$new_seq\n";# print out *.converted.fa
		print O_C ">$new_id\n$new_seq\n"; # print out *.converted.fa
		print O_S ">$new_id\t$new_seq\n";# print out *.converted.SingleLine
	}
	
	close C;
	close A;
}


#
