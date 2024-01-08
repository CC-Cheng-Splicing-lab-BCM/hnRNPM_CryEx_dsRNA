#!/bin/bash

while read line;
do
	left_count=`echo $(echo $line|awk '{print $2}')`
	right_count=`echo $(echo $line|awk '{print $3}')`
	samtools view -b ./CryEx_2023/HTE1_chr21.bam "chr21:$left_count-$right_count" > INTRON
	intron_count=`echo $(samtools view INTRON|awk '(($4*1)>='"$left_count"' && ($4*1) <='"$right_count"')'|wc -l)`
	ex_count=`echo $(echo $line|awk '{print $4}')`
	intron_norm=`echo "scale=2; $intron_count/($right_count - $left_count+100)"|bc`
	len_intron=`echo $(echo $line|awk '{print $3-$2+1}')`
	ex_norm=`echo "scale=2; $ex_count/99"|bc`
	PSI=`echo "scale=2; $intron_norm/($intron_norm + $ex_norm)"|bc`
	printf "$line\t$intron_count\t$len_intron\t$PSI\n" >> chr21_intron.PSI
done < ./CryEx_2023/IR_HTE/SS_count.txt

