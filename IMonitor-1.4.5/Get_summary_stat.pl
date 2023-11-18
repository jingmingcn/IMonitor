#!/usr/bin/perl -w
use strict;

die "perl $0 <dir> <sample> <fa-flag> <sequence-error-flag><CDR3_byconserve>\n" unless(@ARGV==5);

print "#Title\tSeq_num\tRate_of_input(%)\tRate_of_rawdata(%)\n";

my $flag = $ARGV[2];
my $error_f = $ARGV[3];
my $cdr3_bycons_f = $ARGV[4];

my $raw = 0;
if($flag)
{
	open I, "$ARGV[0]/${ARGV[1]}_raw_data_filter_1.txt" or die;
	while(<I>)
	{
	chomp;
	my @line = split;
	if(/All_read_number\(PE=1\):/){
		$line[0] =~s/:$//;
		print "Raw_seq_number(PE=1)\t$line[1]\t-\t-\n";
		$raw = $line[1];
	}
	elsif(/retained_read_num:/){
		my $r_i = sprintf("%0.2f",$line[2]);
		print "Clean_data\t$line[1]\t$r_i\t$r_i\n";
	}
	}
	close I;

	open I, "$ARGV[0]/${ARGV[1]}_merge_filter_1.txt" or die;
	while(<I>)
	{
	chomp;
	my @line = split;
	if(/Merged_seq_raw:/){
		my $r_i = sprintf("%0.2f",$line[2]);
		my $r_r = sprintf("%0.2f",$line[1]/$raw*100);
		print "PE_read_merged\t$line[1]\t$r_i\t",$r_r,"\n";
	}
	if(/Merged_seq_with_high_quality:/){
		my $r_i = sprintf("%0.2f",$line[2]);
		my $r_r = sprintf("%0.2f",$line[1]/$raw*100);
		print "Merged_with_highquality\t$line[1]\t$r_i\t$r_r\n";
	}
	}
	close I;

}


my $before_eff = 0;
open I, "$ARGV[0]/${ARGV[1]}_VDJ.alignment.stat_2.txt" or die;
while(<I>)
{
	chomp;
	my @line = split;

	if($flag==0){
		if(/^sequence_num:/){
			print "Raw_seq_number\t$line[1]\t-\t-\n";
			print "Clean_data\t-\t-\t-\n";
			print "PE_read_merged\t-\t-\t-\n";
			print "Merged_with_highquality\t-\t-\t-\n";
			$raw = $line[1];
		}
	}

	my $r_i = sprintf("%0.2f",$line[2]) if($.>1);
	my $r_r = sprintf("%0.2f",$line[1]/$raw*100);
	if(/^V_alignment_rate:/){
		print "V_alignment\t$line[1]\t$r_i\t$r_r\n";
	}
	elsif(/^D_alignment_rate:/){
		print "D_alignment\t$line[1]\t$r_i\t$r_r\n";
        }
	elsif(/^J_alignment_rate:/){
		print "J_alignment\t$line[1]\t$r_i\t$r_r\n";
	}
	elsif(/^VJ_alignment_rate:/){
		print "VJ_alignment\t$line[1]\t$r_i\t$r_r\n";
		$before_eff = $line[1];
	}


}
close I;
#if(-e "$ARGV[0]/${ARGV[1]}_pcr_seq_error_cor_3.txt" and -s "$ARGV[0]/${ARGV[1]}_pcr_seq_error_cor_3.txt")
if($error_f)
{
	if(-e "$ARGV[0]/${ARGV[1]}_pcr_seq_error_cor_3.txt" and -s "$ARGV[0]/${ARGV[1]}_pcr_seq_error_cor_3.txt")
	{
	open I, "$ARGV[0]/${ARGV[1]}_pcr_seq_error_cor_3.txt" or die;
	my $input;
	while(<I>)
	{
		chomp;
		my @line = split;
		if(/^all_raw_seq:/){
			$input = $line[1];
		}
		if(/new_core_sequence:/){
			my $r_i = sprintf("%0.2f",$line[1]/$input*100);
			my $r_r = sprintf("%0.2f",$line[1]/$raw*100);
			print "PCR_Sequencing_correct\t$line[1]\t$r_i\t$r_r\n";
			$before_eff = $line[1];
		}

	}
	close I;
	}else{
		print "Error: ${ARGV[1]}_pcr_seq_error_cor_3.txt is not existed!\n";
	}
}

open I , "$ARGV[0]/${ARGV[1]}_structure.stat_3.txt" or die;
my $input = $before_eff ;
my $cdr3_found_bycons = 0;


my $stat_raw;

my $flag_cons = 0;
my $eff = 0;
while(<I>)
{
	chomp;
	my @line = split;
	if($cdr3_bycons_f)
	{
		if(/Retained_seq:/){
			$stat_raw = $line[1];
			my $r_i = sprintf("%0.2f",$line[1]/$input*100);
			my $r_r = sprintf("%0.2f",$line[1]/$raw*100);
			print "CDR3_found_VJ\t$line[1]\t$r_i\t$r_r\n";
			$eff = $line[1];
		}
		if(/Find_CDR3_byconserve:/){
			my $r_i = sprintf("%0.2f",$line[1]/$input*100);
			my $r_r = sprintf("%0.2f",$line[1]/$raw*100);
			print "CDR3_found_byconserve\t$line[1]\t$r_i\t$r_r\n";
			$flag_cons = 1;
			$eff = $stat_raw + $line[1];
		}
	}
	else{
		if(/Retained_seq:/){
			my $r_i = sprintf("%0.2f",$line[1]/$input*100);
			my $r_r = sprintf("%0.2f",$line[1]/$raw*100);
			print "CDR3_found_VJ\t$line[1]\t$r_i\t$r_r\n";
			print "CDR3_found_byconserve\t-\t-\t-\n";
			$flag_cons = 1;
			$eff = $line[1];
		}
	}
}
close I;
if($flag_cons == 0){
	print "CDR3_found_byconserve\t0\t0\t0\n";
}
my $r_i2 = sprintf("%0.2f",$eff/$input*100);
my $r_r2 = sprintf("%0.2f",$eff/$raw*100);
print "Effective_data\t$eff\t$r_i2\t$r_r2\n";

print "\n\n";
print "----------Note:--------------\n";
print "Clean_data: filter the Adapter pollution,low quality sequence\n";
print "Effective_data: filter the sequence: 1. cannot find CDR3; 2. V and J strand conflict; 3. CDR3 less than 0bp; 4. sequence abundance filter.\n";


