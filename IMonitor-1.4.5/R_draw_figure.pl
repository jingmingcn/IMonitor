#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin $Script);


=head1
	perl R_draw_figure.pl

=head1 Parameters
	
	-n	sample name
	-i	input files' directory
	-o	out put directory
	-f	the insert-size file
	-m	mutation
	-c	logical value, used for V/J nucleotide composition
	-s	logical value, used for rarefraction curve(Saturation)
	-R	the path for "Rscript" [/opt/blc/genome/bin/Rscript]
	
=cut

my ($name, $compo, $rare_f , $mut_f , $out_dir, $in_dir,$insertsize , $Rscript);

GetOptions(
		"n:s" => \$name,
		"i:s" => \$in_dir,
		"o:s" => \$out_dir,
		"f:s" => \$insertsize,
		"m" => \$mut_f,
		"c" => \$compo,
		"s" => \$rare_f,
		"R:s" => \$Rscript
	  );
die `pod2text $0` if (!$name && !$out_dir && !$in_dir);

$Rscript="/opt/blc/genome/bin/Rscript" if(!defined($Rscript));

my $Bin_p = $Bin;
$Bin = $Bin."/R";

#---------  input files	--------------

`gzip -dc $in_dir/${name}_CDR3_AA.frequency.gz >$in_dir/${name}_CDR3_AA.frequency`;
my $v_usage = "$in_dir/${name}_V_gene.usage";
my $j_usage = "$in_dir/${name}_J_gene.usage";
my $vj_pairing = "$in_dir/${name}_VJ_pairing.usage";
my $cdr3_len = "$in_dir/${name}_CDR3_NT.length";
my $vdj_del = "$in_dir/${name}_VDJ_deletion_nt_len.dis";
my $vdj_ins = "$in_dir/${name}_VDJ_insertion_nt_len.dis";
#my $cdr3_window = "$in_dir/${name}_CDR3_AA_window.stat";
my $cdr3_sec = "$in_dir/${name}_CDR3_AA_section.stat";
my $cdr3_f = "$in_dir/${name}_CDR3_AA.frequency";
my $v_com = "$in_dir/${name}_V_NT_composition.txt";
my $j_com = "$in_dir/${name}_J_NT_composition.txt";
my $rarefract = "$in_dir/${name}_rarefraction_curve.txt";
my $vdj_len = "$in_dir/${name}_VDJ_length_inCDR3.dis";
my $mut_file = "$in_dir/${name}_hypermutation.stat";


#------------
open O, ">$out_dir/tmp" or die;
print O "0\n";
close O;

if(defined($insertsize)){
`$Rscript $Bin/Insert_size_dis.R $insertsize $out_dir/${name}_Insert_size_dis.pdf $name`;
}else{
	$insertsize = "$out_dir/tmp";
}
`$Rscript $Bin/CDR3_length_dis.R $cdr3_len $out_dir/${name}_CDR3_nt_length_dis.pdf $name`;
`$Rscript $Bin/VJ_usage.R $v_usage $out_dir/${name}_V_usage.pdf $name`;
`$Rscript $Bin/VJ_usage.R $j_usage $out_dir/${name}_J_usage.pdf $name`;
`$Rscript $Bin/VDJ_del_length_dis.R $vdj_del $out_dir/${name}_vdj_del_len.pdf $name`;
`$Rscript $Bin/VDJ_ins_length_dis.R $vdj_ins $out_dir/${name}_vdj_ins_len.pdf $name`;
`$Rscript $Bin/VDJ_length_inCDR3_dis.R $vdj_len $out_dir/${name}_vdj_len_inCDR3.pdf $name`;
`$Rscript $Bin/CDR3_freq_dis.R $cdr3_f $cdr3_sec $cdr3_f $out_dir/${name}_CDR3_freq.pdf $name`;
`$Rscript $Bin/VJ_pairing_3d.R $vj_pairing $out_dir/${name}_VJ_pairing_3d $name`;


if($mut_f)
{
	`$Rscript $Bin/Mutation_dis.R $mut_file $out_dir/${name}_mutation_dis.pdf $name`;
}else{
	$mut_file = "$out_dir/tmp";
}


if($rare_f)
{
	`$Rscript $Bin/Rarefraction_curve.R $rarefract $out_dir/${name}_rarefract_curve.pdf $name`;
}else{
	$rarefract = "$out_dir/tmp";
}


if($compo)
{
	my $max_v = `head -1 $v_com`;
	$max_v = (split /\s+/,$max_v)[0];
	my $max_j = `head -1 $j_com`;
	$max_j = (split /\s+/,$max_j)[0];
	`perl $Bin_p/creat_fa_for_weblogo.pl $v_com |$Bin_p/weblogo/seqlogo -f - -F PDF -w 30 -c -k 1 -s $max_v -n -M -t \"Nucleotide Compositon of V($name)\" -o $out_dir/${name}_V_NT_composition`;
	`perl $Bin_p/creat_fa_for_weblogo.pl $j_com |$Bin_p/weblogo/seqlogo -f - -F PDF -w 30 -c -k 1 -s $max_j -n -M -t \"Nucleotide Compositon of J($name)\" -o $out_dir/${name}_J_NT_composition`;
}

`$Rscript $Bin/Overall_plot.R $out_dir/${name}_overall_plot.pdf $name $insertsize $rarefract $cdr3_len $cdr3_f $cdr3_sec $cdr3_f $vdj_len $vdj_del $vdj_ins $mut_file $j_usage $v_usage $vj_pairing`;

unlink("$out_dir/tmp");
`rm $in_dir/${name}_CDR3_AA.frequency`;
