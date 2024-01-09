#!/bin/bash


### This script generates individual subfolder for each expressed exon identified in the sample.
### In this subfolder it contains: 
### junction chromosome coordinates and number of junction reads aligned to that location for both 5'ss and 3'ss of the exon:
### final_junc3SS_HTE_count.txt
### final_junc5SS_HTE_count.txt
cat chr21_exons.coord | while read LINE
do
	echo $LINE
	FILENAME=$LINE
	echo $FILENAME
	FILE=`echo $FILENAME|awk '{print $1}'`
	echo $FILE
	DIRNAME=`echo $FILENAME|awk '{print $1"-"$2}'|sed -r 's/:/-/g'`
	echo $DIRNAME
	DOWN=`echo $FILENAME|awk '{print $1}'|sed 's/.*-//'`
	echo $DOWN
	UP=`echo $FILENAME|awk '{print $1}'|sed 's/.*://'|cut -f1 -d'-'`
	echo $UP
	mkdir ./sub/$DIRNAME
	sed "s/AXX/${FILE}/g" ./CryEx_2023/Code/HTE1/CryEx_Coor_Exact.sh > sub/CryEx.${FILE}.sh
	sed -i "" "s/OXX/${DOWN}/g" sub/CryEx.${FILE}.sh
	sed -i "" "s/UXX/${UP}/g" sub/CryEx.${FILE}.sh
	sed -i "" "s@DXX@/./CryEx_2023/Code/HTE1/sub/${DIRNAME}@g" sub/CryEx.${FILE}.sh
	bash sub/CryEx.${FILE}.sh
	
done
