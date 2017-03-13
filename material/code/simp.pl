#!/usr/bin/perl -w
use strict;

my ($input_file) = $ARGV[0];
my ($output_file1) = $ARGV[1];
my ($output_file2) = $ARGV[2];

my $usage = "This script is to convert summary file into simple file.
summary file:Chromosome Alignment start Alignment end Strand identity Chromosome Alignment start Alignment end score
usage: $0 <summary_file> <simple_file>
";
die $usage if $#ARGV<2;

open(FILE1,$input_file)||die("open $input_file error!\n");
open(FILE2, ">", $output_file1);
open(FILE3, ">", $output_file2);

my @sent=();

while(<FILE1>){
	chomp;
	@sent=split(/\t/,$_);
	if($sent[4]>0.6){
		print FILE2 $sent[0],"\t",$sent[1],"\t",$sent[2],"\t",$sent[3],"\n";
		print FILE3 $_,"\n";
	}

}


close FILE1;
close FILE2;
close FILE3;
exit;
