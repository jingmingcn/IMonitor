#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

unless (@ARGV == 2)
{
	print "Usage:\n";
	print "perl $0 <*.fa>  <out>\n";
	print "";
	exit;
}
my ($in_file1 , $out) = @ARGV;
open IN1 , "$in_file1" or die;
#open IN2 , "$in_file2" or die;
open OUT , ">$out" or die;

my $tolerant_num = 3;

#-----------	read *.fa	-----------------------
my %same_type_temp;
while(<IN1>)
{
	chomp;
	my $id = $_;
	chomp(my $seq = <IN1>);
	if((split /:/ , $id)[-1] eq "p") # "p" print out directly
	{
		print OUT "$id\n$seq\n";
	}
	else
	{
		my $type = (split /:/,$id)[0];
		$same_type_temp{$type}{(split /:/ , $id)[3]} = "$id\n$seq\n";
	}
}
close IN1;

# conduct the "t" sequences
for my $type(keys %same_type_temp)
{
	my %combined_inf;
	my %seq_backup;

	# use 2 "for" to combined sequences  
	for my $len(sort {$b<=>$a} keys %{$same_type_temp{$type}})
	{
		next if(exists $seq_backup{$len}); # 

		my $primer_len = (split /:/ , $same_type_temp{$type}{$len})[4];

		for my $len2 (sort {$b<=>$a} keys %{$same_type_temp{$type}})	
		{
			next if($len <= $len2); # 
			if($len-$len2 <= 3)# combine to one sequence during less than 3bp length difference
			{
				my $primer_len2 = (split /:/ , $same_type_temp{$type}{$len2})[4]+$len-$len2; # primer len + sequences length difference
				$primer_len = $primer_len2 if($primer_len < $primer_len2);
				$seq_backup{$len2} = 1;		
			}
		}
		my @temp = split /:/, $same_type_temp{$type}{$len};
		$temp[4] = $primer_len; # change the primer length
	
		$temp[0] = join ":" , @temp;
		print OUT "$temp[0]";
		$seq_backup{$len} = 1;
	}
}
