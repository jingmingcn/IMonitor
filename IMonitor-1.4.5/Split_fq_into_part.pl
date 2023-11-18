#!/usr/bin/perl -w 
use strict;

die "perl $0 <in.gz> <> <dir> <sample><num>\n" unless(@ARGV==5);

my $num = $ARGV[4];
#my $all = `awk '\$1=="erged_seq_with_high_quality:"{print}' $ARGV[1]`;print "$ARGV[1]\n$all\n";
open I, "$ARGV[1]" or die;
my $all;
while(<I>)
{
	if(/Unique_seq_num:/){
		chomp;
		$all = (split)[1];
	}
}
close I;

my $interval = int $all/$num;

if($ARGV[0]=~/\.gz$/){
	open I, "gzip -dc $ARGV[0]|" or die;
}else{
	open I, "$ARGV[0]" or die;
}


my $flag = 0;
my $all_f = 0;
my %new;
my $last = $interval*($num-1);

my %back;

while(<I>)
{
	chomp;
	s/^>//;
	my $id = $_;
	my $raw = $id;
	$id =~ s/^cp//;
	my $abund;
	($id,$abund) = split /c:/,$id;
	$back{$id} = $raw;
	chomp(my $seq = <I>);
	$flag++;
	$all_f++;
	$new{$id} = $seq;

	if($flag==$interval && $all_f<=$last)
	{
		my $n = int $all_f/$interval;
		open O, ">$ARGV[2]/$ARGV[3].fa.tmp.$n" or die;
		for(sort {$a<=>$b} keys %new)
		{
			print O ">$back{$_}\n$new{$_}\n";
		}
		close O;
		%new = ();
		%back = ();
		$flag = 0;
	}
}
close I;


open O, ">$ARGV[2]/$ARGV[3].fa.tmp.$num" or die;
for(sort {$a<=>$b} keys %new)
{
	print O ">$back{$_}\n$new{$_}\n";
}
close O;
