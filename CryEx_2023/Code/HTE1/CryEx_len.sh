#!/bin/bash

### This script generates the number of splice sites identified for both 5'SS and 3'SS in any expressed exons in a sample 
### (preferablly 1 and 1 for 5'ss and 3'ss separately, which is the majority luckily. There are situations when there is no 5'ss or no 3'ss or multiple ss)
touch CryEx_coord_len.bed
cd ./CryEx_2023/Code/HTE1/sub/
for d in */
do
	echo $d
	line=$d
	### reading into the splicing/junction information for both ss in an expressed exon in a sample
	SS3_file=`printf "./CryEx_2023/Code/HTE1/sub/${line}final_junc3SS_HTE_count.txt"`
	SS5_file=`printf "./CryEx_2023/Code/HTE1/sub/${line}final_junc5SS_HTE_count.txt"`
	chr=`echo $line|cut -d'-' -f 1`
	### counting how many junction sites there representing 3'SS
	SS3=`echo $(cat $SS3_file|wc -l)`
	echo $SS3
	### counting how many junction sites there representing 3'SS
	SS5=`echo $(cat $SS5_file|wc -l)`
	echo $SS5
	
	
	printf "$line\t$chr\t$SS3\t$SS5\n" >> ../CryEx_coord_len.bed

done
exit	

