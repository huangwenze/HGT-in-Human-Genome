#!/usr/bin/perl -w
use strict;

my ($input_file) = $ARGV[0];
my ($output_file) = $ARGV[1];

my $usage = "This script is to convert axt format file into summary file and record the distribution of identity.
summary file:Chromosome Alignment start Alignment end Strand identity Chromosome Alignment start Alignment end score
usage: $0 <axt_format_file> <summary_file>
";

die $usage if $#ARGV<1;

open(FILE1,$input_file)||die("open $input_file error!\n");
open(FILE2, ">", $output_file);
open(FILE3, ">>count.txt");

my $sen1=" ";
my $sen2=" ";
my @sent=();
my @sent1=();
my @sent2=();
my $i=0;
my $j=0;
my $same_num=0;
my $perc=0.0;
my $r=0;
my $flag=0;
my @dist=();

for($j=0; $j<=9; $j++){
    @dist=(@dist,0);
}
 

while(<FILE1>){
	chomp;
	if($flag==3){
	for($i=0; $i<=$#sent1 && $i<=$#sent2 ; $i++){
		if($sent1[$i] eq $sent2[$i]){
			$same_num++;
		}
	}
	$perc=$same_num/($#sent1+1);
	for($r=0;$r<=9;$r++){
	    if($perc<=$r*0.1+0.1){
		$dist[$r]++;
		last;
	    }
	}
	print FILE2 $sent[1],"\t",$sent[2],"\t",$sent[3],"\t",$sent[7],"\t", sprintf"%.3f\t",$perc;
	print FILE2 $sent[4],"\t",$sent[5],"\t",$sent[6],"\t",$sent[8],"\n";
	$flag=0;
	$same_num=0;
	$perc=0.0;
	}
	if(/^[0-9]+\schr[0-9XY]+\s[0-9]+\s[0-9]+\s/){
	    @sent=split(/\s/,$_);
 
	    @sent1=();
	    @sent2=();
	    $flag=1;
	}
	elsif($flag==1){
	    $_=uc($_);
	    @sent1=split(//,$_);
	    $flag=2;
	}
	elsif($flag==2){
	    $_=uc($_);
	    @sent2=split(//,$_);
	    $flag=3;
	}
	
}

if($flag==3){
	for($i=0; $i<=$#sent1 && $i<=$#sent2 ; $i++){
		if($sent1[$i] eq $sent2[$i]){
			$same_num++;
		}
	}
	$perc=$same_num/($#sent1+1);
	for($r=0;$r<=9;$r++){
	    if($perc<=$r*0.1+0.1){
		$dist[$r]++;
		last;
	    }
	}
	print FILE2 $sent[1],"\t",$sent[2],"\t",$sent[3],"\t",$sent[7],"\t", sprintf"%.3f\t",$perc;
	print FILE2 $sent[4],"\t",$sent[5],"\t",$sent[6],"\t",$sent[8],"\n";
	$flag=0;
	$same_num=0;
	$perc=0.0;
}

print FILE3 "convert ", $input_file, " to ", $output_file, "\n";
print FILE3 "distribution:";
for($j=0; $j<=9; $j++){
    print FILE3 "\t",$dist[$j];
}
print FILE3 "\n";

close FILE1;
close FILE2;
close FILE3;
exit;
