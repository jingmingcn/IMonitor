#!/usr/bin/perl -w 
use strict;

die "perl $0 <in.gz> <*.backup> <out.gz>\n" unless(@ARGV==3);

my %id;
open I, "$ARGV[1]" or die;
while(<I>)
{
	chomp;
	my @line = split;
	if($line[0]=~/:/){
		my @t = split /:/,$line[0];
		$id{$line[1]} = \@t;
	}else{
		$id{$line[1]}->[0] = $line[0];
	}
}
close I;


open I, "gzip -dc $ARGV[0]|" or die;
open O, "|gzip > $ARGV[2]" or die;
my $head = <I>;
print O "$head";
while(<I>)
{
	chomp;
	my @line = split;
	$line[2] = &change_id($line[2]) if($line[2] ne "NA");
	$line[3] = &change_id($line[3]) if($line[3] ne "NA");
	$line[4] = &change_id($line[4]) if($line[4] ne "NA");
	my $new = join "\t" , @line;
	print O "$new\n";
}
close I;

sub change_id
{
	my $old = shift @_;
	my $new;
	if(!exists $id{$old}){
		print "Error: the id $old is not exist!\n";
		exit;
	}

	if(scalar @{$id{$old}} ==1 ){
		return $id{$old}->[0];
	}
	else # multiple sequences are the same, then select it at random
	{
		my $rand = int rand(scalar @{$id{$old}});
		return $id{$old}->[$rand];
	}
}

