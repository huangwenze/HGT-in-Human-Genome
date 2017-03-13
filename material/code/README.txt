axttofile.pl : perl script to convert axt format file into summary file and record the distribution of identity.
Usage is shown as follows : 
perl axttofile.pl <axt_format_file> <summary_file> 
eg: perl axttofile.pl chr1.hg19.panTro4.net.axt Chimp_chr1_m

simp.pl : perl script to convert summary file into simple file.
Usage is shown as follows :
perl simp.pl <summary_file> <simple_file>
eg: perl simp.pl Chimp_chr8_m Chimp_chr8_s Chimp_chr8_l

prepare.pl : perl script to generate a record_file to record the comparison information.
Usage is shown as follows :
perl prepare.pl

select6 : C script to compare two DNA fragment position and record the overlapping fragment in record_file.txt
Usage is shown as follows:
./select6 -c -a start,end -b start,end <gene_position_file1> <gene_position_file2>
eg: ./select6 -c -a 2,3 -b 2,3 Zebra_finch_chr1_s Chimp_chr1_s
mv record_file.txt Zebra_finch_chr1_f

result.pl : perl script to get the information from comparison file and generate result in Zebra_finch_chr1_r.
Usage is shown as follows:
perl result.pl <input_file> <input_data> <output_file>
eg: perl result.pl Zebra_finch_chr1_f Zebra_finch_chr1_l Zebra_finch_chr1_r
