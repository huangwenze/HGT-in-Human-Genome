#!/usr/bin/perl -w
use strict;

my ($input_file) = $ARGV[0];
my ($input_data) = $ARGV[1];
my ($output_file) = $ARGV[2];
my $usage = "This script is to get the information from comparison file.
usage: $0 <input_file> <input_data> <output_file>
";
die $usage if $#ARGV<2;

open(INPFILE,$input_file)||die("open $input_file error!\n");
open(INPDATA,$input_data)||die("open $input_data error!\n");
open(FILE, ">", $output_file);
my @info=();
my @numb=();
my $sen1="";
my $sen2="";
my $sen3="";
my $i=0;
my $j=0;
my $r=0;
my $flag=0;
while($sen2=<INPDATA>){
	$sen1=<INPFILE>;
    chomp($sen1);
	chomp($sen2);
	$j=0;
	$r=0;
	@info = split(//,$sen1);
	@numb = split(/\t/,$sen2);
	for($i=0;$i<=41;$i++){
		if($info[$i] eq 'n'){
			$j++;
		}
	}
    for($i=$i;$i<=$#info;$i++){                                                
	    if($info[$i] eq 'n'){
		$r++;
	    }
    }
	if($j <= 8  && $r > 0 && $numb[2]-$numb[1] >= 1000 ){		
    print FILE $sen2,"\t",$j,"\t",$r,"\t",$sen1,"\n";
	}
}

close INPFILE;
close INPDATA;
close FILE;
exit;
