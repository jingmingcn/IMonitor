#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

unless (@ARGV == 3)
{
	print "Usage:\n";
	print "perl $0 <V/J.fa><start_postion><out>\n";
	print "";
	exit;
}
my ($in_file1 , $start_pos , $out) = @ARGV;
open IN1 , "$in_file1" or die;
#open IN2 , "$in_file2" or die;
open OUT , ">$out" or die;


my %seq_new;
my $max_len = 0;
my $numbering = 0;
$/=">";
<IN1>;
while(<IN1>)
{
	chomp;
	my @line = split /\n+/,$_;
	my $id = (split /\|/,$line[0])[1];
	if((split /\|/, $line[0])[3] ne "F" and (split /\|/, $line[0])[3] ne "(F)" and (split /\|/, $line[0])[3] ne "[F]"){$id.="_unF";}# for pseudogene , add a tag
	shift @line;
	my $raw_seq="@line";
	$raw_seq=~s/\s+//g;

	$raw_seq =~ tr/acgt./ACGT-/;

	next if(length($raw_seq)<$start_pos);
	my $sub_seq = substr($raw_seq , $start_pos-1);# cut sequence
	$max_len = length($sub_seq) if($max_len < length($sub_seq));

	my $sub_seq_only = $sub_seq; # get bases length
	$sub_seq_only =~ s/-//g; 
	my $len = length($sub_seq_only);

#	print OUT ">$id:$numbering:$start_pos:$len:0:0:t\n$sub_seq\n";
	$id = ">$id:$numbering:0:$len:0:0:0:t";
	$seq_new{$id} = $sub_seq;
	$numbering++;
}
close IN1;
$/="\n";

for (keys %seq_new)
{
	my $len = $max_len - length($seq_new{$_});
	my $seq = $seq_new{$_}.("-" x $len);
	print OUT "$_\n$seq\n";
}

