#!/usr/bin/perl -w
use strict;

my $usage = "This script is to generate a record_file to record the comparison information.
usage: $0 
";
die $usage if $#ARGV<-1;

open(FILE, ">record_file.txt");
my $i = 0;
for($i=0; $i<10000000; $i++){
	print FILE ">\n";
}
close FILE;
exit;
