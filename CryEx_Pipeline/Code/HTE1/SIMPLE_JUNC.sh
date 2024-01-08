#!/bin/bash

### This script calculates for PSI values for all combinations of 5'ss and 3'ss found in any expressed exon identified in a sample: 
### In the end, each well-defined exon with its PSI value quantification is categorized into three types:
### FIRST: with only 5'ss but no 3'ss (FIRST.PSI)
### MIDDLE: having both 5'ss and 3'ss (MIDDLE.PSI)
### LAST: with only 3'ss but no 5'ss (LAST.PSI)

cat ./CryEx_2023/HTE1_chr21_transcripts.gtf|grep "exon"|awk '($1=="chr21")' > ./CryEx_2023/Code/HTE1/chr21_exons.gtf
cat ./CryEx_2023/Code/HTE1/chr21_exons.gtf|awk '{print $1"-"$4"-"$5"-"$12"-"$14"/""\t"$7}'|tr -d '"'|tr -d ';' > strandedness
sort -k1,1 strandedness > strandedness.sort
sort -k1,1 CryEx_coord_len.bed > CryEx_coord_len.bed.sort
#join -j1 --nocheck-order CryEx_coord_len.bed.sort strandedness.sort > CryEx_coord_len.bed.strand
# For MacOS join command run:
join CryEx_coord_len.bed.sort strandedness.sort > CryEx_coord_len.bed.strand
perl -p -i -e 's/ /\t/g' CryEx_coord_len.bed.strand
cat CryEx_coord_len.bed.strand |awk '($3!=0 && $4!=0)'> middle.coord

cat CryEx_coord_len.bed.strand |awk '($3!=0 && $4==0)'|awk '($5=="-")'> first_1.coord
cat CryEx_coord_len.bed.strand |awk '($3==0 && $4!=0)'|awk '($5=="+")'> first_2.coord
cat first_1.coord first_2.coord > first.coord

cat CryEx_coord_len.bed.strand |awk '($3==0 && $4!=0)'|awk '($5=="-")'> last_1.coord
cat CryEx_coord_len.bed.strand |awk '($3!=0 && $4==0)'|awk '($5=="+")'> last_2.coord
cat last_1.coord last_2.coord > last.coord

cat CryEx_coord_len.bed.strand |awk '($3==0 && $4==0)'> notspliced.coord

rm strandedness
rm strandedness.sort
rm CryEx_coord_len.bed.sort
rm CryEx_coord_len.bed.strand
rm first_1.coord
rm first_2.coord
rm last_1.coord
rm last_2.coord





### FIRST.bed
cat first.coord |awk '{print $1}'|cut -f1 -d '-' > chr
cat first.coord |awk '{print $1}'|cut -f2 -d '-' > start
cat first.coord |awk '{print $1}'|cut -f3 -d '-' > end
cat first.coord |awk '{print $5}'>sign
paste -d ':' chr start end sign > lable
paste -d '\t' lable first.coord|awk '!seen[$1]++'|awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$5"\t"$8}' > trim_first.coord

while read line;
do
	echo $line
	sign=`echo $(echo $line|awk '{print $5}')`
	path=`echo $(echo $line|awk '{print $1}')`
	chr=`echo $(echo $line|awk '{print $2}')`
	
	
	if [[ "$sign" == "+" ]]
	then 
		cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc5SS_HTE_count.txt > SS5_file	
		SS3=`echo $path|cut -d'-' -f 2`	
		while read ss;
		do 
			echo $ss
			SS5=`echo $ss|awk '{print $1}'`
			cat add_SS_count.txt|grep $SS5|awk '($3> '"$SS5"')'> DN
			while read ll;
			do
				echo $ll
				DN_coord=`echo $ll|awk '{print $3}'`
				DN_count=`echo $ll|awk '{print $4}'`
			
				cat add_SS_count.txt|awk '($3 == '"$DN_coord"')'|awk '($2< '"$SS3"')'> SK
				tmp=`echo $(cat SK)`
				if [[ ! -z "$tmp" ]]
				then
					
					while read kk;
					do
						echo $kk
						SK_UP=`echo $kk|awk '{print $2}'`
					
						SK_UP_count=`echo $kk|awk '{print $4}'`
						SK_UP_coord=$SK_UP
		
					done < SK
				elif [ -z "$tmp" ]
				then 
					SK_UP_count=0
					SK_UP_coord=0
			

				else
					echo "what?"
				fi
				samtools view -b /Users/rongzheng/Downloads/CryEx_2023/HTE1_chr21_sorted.bam "chr21:$SS3-$SS5"> GENE	
				Ex_count=`echo $(samtools view GENE|awk '(($4*1)>='"$SS3"' && ($4*1)<='"$SS5"')'|awk '($6=="100M")'|wc -l)`
				add_2=`echo "scale=10; ($DN_count+$Ex_count)/($SS5-$SS3+100)"|bc`
				tot_2=`echo "scale=10; $add_2+($SK_UP_count/99)"|bc`
				PSI=`echo "scale=2; $add_2/$tot_2"|bc`
			printf "$chr\t$SS3\t$SS5\t$SK_UP_coord\t$DN_coord\t$SK_UP_count\t$DN_count\t$Ex_count\t$PSI\t$sign\n" >> first_plus.PSI
			
			done < DN
		done < SS5_file
		
		
	elif [ "$sign" == "-" ]
	then
		cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc3SS_HTE_count.txt > SS3_file
		SS5=`echo $path|cut -d'-' -f 3`
		while read ss;
		do
			echo $ss
			SS3=`echo $ss|awk '{print $1}'`
			cat add_SS_count.txt|grep $SS3|awk '($2< '"$SS3"')'> UP
			while read ll;
			do
				echo $ll
				
				UP_coord=`echo $ll|awk '{print $2}'`
				UP_count=`echo $ll|awk '{print $4}'`
			
				cat add_SS_count.txt|awk '($2 == '"$UP_coord"')'|awk '($3> '"$SS5"')'> SK
				tmp=`echo $(cat SK)`
				if [[ ! -z "$tmp" ]]
				then
					while read kk;
					do 
						echo $kk
						SK_DN=`echo $kk|awk '{print $3}'`
						SK_DN_count=`echo $kk|awk '{print $4}'`
						SK_DN_coord=$SK_DN
					done < SK
			
				elif [ -z "$tmp" ]
				then 
					SK_DN_count=0
					SK_DN_coord=0
						
			

				else
					echo "what?"
				fi
		
				samtools view -b /Users/rongzheng/Downloads/CryEx_2023/HTE1_chr21_sorted.bam "chr21:$SS3-$SS5"> GENE	
				Ex_count=`echo $(samtools view GENE|awk '(($4*1)>='"$SS3"' && ($4*1)<='"$SS5"')'|awk '($6=="100M")'|wc -l)`
				add_2=`echo "scale=10; ($UP_count+$Ex_count)/($SS5-$SS3+100)"|bc`
				tot_2=`echo "scale=10; $add_2+($SK_DN_count/99)"|bc`
				PSI=`echo "scale=2; $add_2/$tot_2"|bc`
			printf "$chr\t$SS3\t$SS5\t$UP_coord\t$SK_DN_coord\t$UP_count\t$SK_DN_count\t$Ex_count\t$PSI\t$sign\n" >> first_minus.PSI
				
			

			done < UP
		done < SS3_file
			
		
		
	else
		echo "what?"
	fi
done < first.coord
cat first_minus.PSI first_plus.PSI|awk '!seen[$0]++' > first.PSI
#rm first_minus.PSI
#rm first_plus.PSI
rm DN
rm SK
rm UP
rm SS3_file
rm SS5_file


### LAST.bed
cat last.coord |awk '{print $1}'|cut -f1 -d '-' > chr
cat last.coord |awk '{print $1}'|cut -f2 -d '-' > start
cat last.coord |awk '{print $1}'|cut -f3 -d '-' > end
cat last.coord |awk '{print $5}'>sign
paste -d ':' chr start end sign > lable
paste -d '\t' lable last.coord|awk '!seen[$1]++'|awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$5"\t"$8}' > trim_last.coord

while read line;
do
	echo $line
	sign=`echo $(echo $line|awk '{print $5}')`
	path=`echo $(echo $line|awk '{print $1}')`
	chr=`echo $(echo $line|awk '{print $2}')`
	
	
	if [[ "$sign" == "+" ]]
	then 
		cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc3SS_HTE_count.txt > SS3_file
		SS5=`echo $path|cut -d'-' -f 3`	
		while read ss;
		do 
			echo $ss
			SS3=`echo $ss|awk '{print $1}'`
			cat add_SS_count.txt|grep $SS3|awk '($2< '"$SS3"')'> UP
			while read ll;
			do
				echo $ll
				UP_coord=`echo $ll|awk '{print $2}'`
				UP_count=`echo $ll|awk '{print $4}'`
			
				cat add_SS_count.txt|awk '($2 == '"$UP_coord"')'|awk '($3> '"$SS5"')'> SK
				tmp=`echo $(cat SK)`
				if [[ ! -z "$tmp" ]]
				then
					
					while read kk;
					do
						echo $kk
						SK_DN=`echo $kk|awk '{print $3}'`
					
						SK_DN_count=`echo $kk|awk '{print $4}'`
						SK_DN_coord=$SK_DN
		
					done < SK
				elif [ -z "$tmp" ]
				then 
					SK_DN_count=0
					SK_DN_coord=0
			

				else
					echo "what?"
				fi
				samtools view -b /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/HTE1_chr21_sorted.bam "chr21:$SS3-$SS5"> GENE	
				Ex_count=`echo $(samtools view GENE|awk '(($4*1)>='"$SS3"' && ($4*1)<='"$SS5"')'|awk '($6=="100M")'|wc -l)`
				add_2=`echo "scale=10; ($UP_count+$Ex_count)/($SS5-$SS3+100)"|bc`
				tot_2=`echo "scale=10; $add_2+($SK_DN_count/99)"|bc`
				PSI=`echo "scale=2; $add_2/$tot_2"|bc`
			printf "$chr\t$SS3\t$SS5\t$UP_coord\t$SK_DN_coord\t$UP_count\t$SK_DN_count\t$Ex_count\t$PSI\t$sign\n" >> last_plus.PSI
			
			done < UP
	
		
		done < SS3_file
	elif [ "$sign" == "-" ]
	then
		cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc5SS_HTE_count.txt > SS5_file
		SS3=`echo $path|cut -d'-' -f 2`
		while read ss;
		do
			echo $ss
			SS5=`echo $ss|awk '{print $1}'`
			cat add_SS_count.txt|grep $SS5|awk '($3> '"$SS5"')'> DN
			while read ll;
			do
				echo $ll
				
				DN_coord=`echo $ll|awk '{print $3}'`
				DN_count=`echo $ll|awk '{print $4}'`
			
				cat add_SS_count.txt|awk '($3 == '"$DN_coord"')'|awk '($2< '"$SS3"')'> SK
				tmp=`echo $(cat SK)`
				if [[ ! -z "$tmp" ]]
				then
					while read kk;
					do 
						echo $kk
						SK_UP=`echo $kk|awk '{print $2}'`
						SK_UP_count=`echo $kk|awk '{print $4}'`
						SK_UP_coord=$SK_UP
					done < SK
			
				elif [ -z "$tmp" ]
				then
					SK_UP_count=0
					SK_UP_coord=0
						
			

				else
					echo "what?"
				fi
		
				samtools view -b /Users/rongzheng/Downloads/CryEx_2023/HTE1_chr21_sorted.bam "chr21:$SS3-$SS5"> GENE	
				Ex_count=`echo $(samtools view GENE|awk '(($4*1)>='"$SS3"' && ($4*1)<='"$SS5"')'|awk '($6=="100M")'|wc -l)`
				add_2=`echo "scale=10; ($DN_count+$Ex_count)/($SS5-$SS3+100)"|bc`
				tot_2=`echo "scale=10; $add_2+($SK_UP_count/99)"|bc`
				PSI=`echo "scale=2; $add_2/$tot_2"|bc`
			printf "$chr\t$SS3\t$SS5\t$SK_UP_coord\t$DN_coord\t$SK_UP_count\t$DN_count\t$Ex_count\t$PSI\t$sign\n" >> last_minus.PSI
				
			

			done < DN
		done < SS5_file
		
	else
		echo "what?"
	fi
done < trim_last.coord
cat last_minus.PSI last_plus.PSI|awk '!seen[$0]++' > last.PSI
#rm last_minus.PSI
#rm last_plus.PSI



### MIDDLE.bed
cat middle.coord |awk '{print $1}'|cut -f1 -d '-' > chr
cat middle.coord |awk '{print $1}'|cut -f2 -d '-' > start
cat middle.coord |awk '{print $1}'|cut -f3 -d '-' > end
cat middle.coord |awk '{print $5}'>sign
paste -d ':' chr start end sign > lable
paste -d '\t' lable middle.coord|awk '!seen[$1]++'|awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6}' > trim_middle.coord

while read line;
do
	echo $line
	sign=`echo $(echo $line|awk '{print $5}')`
	path=`echo $(echo $line|awk '{print $1}')`
	chr=`echo $(echo $line|awk '{print $2}')`
	
	
	if [[ "$sign" == "+" ]]
	then 
		
		min_SS3=`echo $(cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc3SS_HTE_count.txt|sort -k1,1n|head -1)`
		min_SS3_coord=`echo $min_SS3|awk '{print $1}'`
		min_SS5=`echo $(cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc5SS_HTE_count.txt|awk '($1< '"$min_SS3_coord"')')`
		if [[ ! -z "$min_SS5" ]]
		then
			
			SS5=`echo $min_SS5|awk '{print $1}'`
			SS3=`echo $path|cut -d'-' -f 2`
			cat add_SS_count.txt|grep $SS5|awk '($3> '"$SS5"')' > DN
			while read zz;
			do
				echo $zz
				DN_coord=`echo $zz|awk '{print $3}'`
				DN_count=`echo $zz|awk '{print $4}'`
				mvp=`echo $(cat add_SS_count.txt|grep $DN_coord|awk '($2< '"$SS3"')')`
				if [[ ! -z "$mvp" ]]
				then
					echo $yy
					UP_coord=`echo $mvp|awk '{print $2}'`
					UP_count=0
					SK_count=`echo $mvp|awk '{print $4}'`
				elif [ -z "$mvp" ]
				then 
					UP_coord=0
					UP_count=0
					SK_count=0
				else
					echo "what?"
				fi
				samtools view -b /Users/rongzheng/Downloads/CryEx_2023/HTE1_chr21_sorted.bam "chr21:$SS3-$SS5"> GENE	
				Ex_count=`echo $(samtools view GENE|awk '(($4*1)>='"$SS3"' && ($4*1)<='"$SS5"')'|awk '($6=="100M")'|wc -l)`
				add_2=`echo "scale=10; ($UP_count+$Ex_count+$DN_count)/($SS5-$SS3+100)"|bc`
				tot_2=`echo "scale=10; $add_2+($SK_count/99)"|bc`
				PSI=`echo "scale=2; $add_2/$tot_2"|bc`
				printf "$chr\t$SS3\t$SS5\t$UP_coord\t$DN_coord\t$UP_count\t$DN_count\t$SK_count\t$Ex_count\t$PSI\t$sign\t$path\n" >> middle_first_plus.PSI
			done < DN
		elif [ -z "$min_SS5" ]
		then
			
			echo "noSS5only"
		else
			echo "what?"
		fi
		max_SS5=`echo $(cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc5SS_HTE_count.txt|sort -k1,1n|tail -1)`
		max_SS5_coord=`echo $max_SS5|awk '{print $1}'`
		max_SS3=`echo $(cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc3SS_HTE_count.txt|awk '($1> '"$max_SS5_coord"')')`
		if [[ ! -z "$max_SS3" ]]
		then
			
			SS3=`echo $max_SS3|awk '{print $1}'`
			SS5=`echo $path|cut -d'-' -f 3`
			cat add_SS_count.txt|grep $SS3|awk '($2< '"$SS3"')' > UP
			while read zz;
			do
				echo $zz
				UP_coord=`echo $zz|awk '{print $2}'`
				UP_count=`echo $zz|awk '{print $4}'`
				mvp=`echo $(cat add_SS_count.txt|grep $UP_coord|awk '($3> '"$SS5"')')`
				if [[ ! -z "$mvp" ]]
				then
					echo $yy
					DN_coord=`echo $mvp|awk '{print $3}'`
					DN_count=0
					SK_count=`echo $mvp|awk '{print $4}'`
				elif [ -z "$mvp" ]
				then 
					DN_coord=0
					DN_count=0
					SK_count=0
				else
					echo "what?"
				fi
				samtools view -b /Users/rongzheng/Downloads/CryEx_2023/HTE1_chr21_sorted.bam "chr21:$SS3-$SS5"> GENE	
				Ex_count=`echo $(samtools view GENE|awk '(($4*1)>='"$SS3"' && ($4*1)<='"$SS5"')'|awk '($6=="100M")'|wc -l)`
				add_2=`echo "scale=10; ($UP_count+$Ex_count+$DN_count)/($SS5-$SS3+100)"|bc`
				tot_2=`echo "scale=10; $add_2+($SK_count/99)"|bc`
				PSI=`echo "scale=2; $add_2/$tot_2"|bc`
				printf "$chr\t$SS3\t$SS5\t$UP_coord\t$DN_coord\t$UP_count\t$DN_count\t$SK_count\t$Ex_count\t$PSI\t$sign\t$path\n" >> middle_last_plus.PSI
			done < UP
		elif [ -z "$max_SS3" ]
		then
			
			echo "noSS5only"
		else
			echo "what?"
		fi
						
			
		cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc3SS_HTE_count.txt > SS3_file
		
		
		while read ss;
		do	
			echo $ss
			SS3=`echo $ss|awk '{print $1}'`
			cat add_SS_count.txt|grep $SS3|awk '($2< '"$SS3"')' > UP
			while read ll;
			do
				echo $ll
				UP_coord=`echo $ll|awk '{print $2}'`
				UP_count=`echo $ll|awk '{print $4}'`
				cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc5SS_HTE_count.txt|awk '($1> '"$SS3"')' > SS5_file
				sud=`echo $(cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc5SS_HTE_count.txt|awk '($1> '"$SS3"')')`
				if [[ ! -z "$sud" ]]
				then
				
					while read mm;
					do
						echo $mm
						SS5=`echo $mm|awk '{print $1}'`
						cat add_SS_count.txt|grep $SS5|awk '($3> '"$SS5"')'> DN
						while read nn;
						do
							echo $nn
							DN_coord=`echo $nn|awk '{print $3}'`
							DN_count=`echo $nn|awk '{print $4}'`
							tmp=`echo $(cat add_SS_count.txt|awk '($2=='"$UP_coord"')'|awk '($3=='"$DN_coord"')')` 
					
							if [[ ! -z "$tmp" ]]
							then
								echo $tmp
								SK_count=`echo $tmp|awk '{print $4}'`
							elif [ -z "$tmp" ]
							then
								SK_count=0
							else
								echo "what?"
							fi
							samtools view -b /Users/rongzheng/Downloads/CryEx_2023/HTE1_chr21_sorted.bam "chr21:$SS3-$SS5"> GENE	
							Ex_count=`echo $(samtools view GENE|awk '(($4*1)>='"$SS3"' && ($4*1)<='"$SS5"')'|awk '($6=="100M")'|wc -l)`
							add_2=`echo "scale=10; ($UP_count+$Ex_count+$DN_count)/($SS5-$SS3+100)"|bc`
							tot_2=`echo "scale=10; $add_2+($SK_count/99)"|bc`
							PSI=`echo "scale=2; $add_2/$tot_2"|bc`
						printf "$chr\t$SS3\t$SS5\t$UP_coord\t$DN_coord\t$UP_count\t$DN_count\t$SK_count\t$Ex_count\t$PSI\t$sign\t$path\n" >> middle_plus.PSI

						done < DN

					done < SS5_file
				elif [ -z "$sud" ]
				then
					echo "didbefore"
				else
					echo "what?"
				fi
				
			done < UP
			
		done < SS3_file
		
		
	elif [ "$sign" == "-" ]
	then
		
		min_SS3=`echo $(cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc3SS_HTE_count.txt|sort -k1,1n|head -1)`
		min_SS3_coord=`echo $min_SS3|awk '{print $1}'`
		min_SS5=`echo $(cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc5SS_HTE_count.txt|awk '($1<'"$min_SS3_coord"')')`
		if [[ ! -z "$min_SS5" ]]
		then
			
			SS5=`echo $min_SS5|awk '{print $1}'`
			SS3=`echo $path|cut -d'-' -f 2`
			cat add_SS_count.txt|grep $SS5|awk '($3> '"$SS5"')' > DN
			while read zz;
			do
				echo $zz
				DN_coord=`echo $zz|awk '{print $3}'`
				DN_count=`echo $zz|awk '{print $4}'`
				mvp=`echo $(cat add_SS_count.txt|grep $DN_coord|awk '($2< '"$SS3"')')`
				if [[ ! -z "$mvp" ]]
				then
					echo $yy
					UP_coord=`echo $mvp|awk '{print $2}'`
					UP_count=0
					SK_count=`echo $mvp|awk '{print $4}'`
				elif [ -z "$mvp" ]
				then 
					UP_coord=0
					UP_count=0
					SK_count=0
				else
					echo "what?"
				fi
				samtools view -b /Users/rongzheng/Downloads/CryEx_2023/HTE1_chr21_sorted.bam "chr21:$SS3-$SS5"> GENE	
				Ex_count=`echo $(samtools view GENE|awk '(($4*1)>='"$SS3"' && ($4*1)<='"$SS5"')'|awk '($6=="100M")'|wc -l)`
				add_2=`echo "scale=10; ($UP_count+$Ex_count+$DN_count)/($SS5-$SS3+100)"|bc`
				tot_2=`echo "scale=10; $add_2+($SK_count/99)"|bc`
				PSI=`echo "scale=2; $add_2/$tot_2"|bc`
				printf "$chr\t$SS3\t$SS5\t$UP_coord\t$DN_coord\t$UP_count\t$DN_count\t$SK_count\t$Ex_count\t$PSI\t$sign\t$path\n" >> middle_last_minus.PSI
			done < DN
		elif [ -z "$min_SS5" ]
		then
			
			echo "noSS5only"
		else
			echo "what?"
		fi
		max_SS5=`echo $(cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc5SS_HTE_count.txt|sort -k1,1n|tail -1)`
		max_SS5_coord=`echo $max_SS5|awk '{print $1}'`
		max_SS3=`echo $(cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc3SS_HTE_count.txt|awk '($1> '"$max_SS5_coord"')')`
		if [[ ! -z "$max_SS3" ]]
		then
			
			SS3=`echo $max_SS3|awk '{print $1}'`
			SS5=`echo $path|cut -d'-' -f 3`
			cat add_SS_count.txt|grep $SS3|awk '($2< '"$SS3"')' > UP
			while read zz;
			do
				echo $zz
				UP_coord=`echo $zz|awk '{print $2}'`
				UP_count=`echo $zz|awk '{print $4}'`
				mvp=`echo $(cat add_SS_count.txt|grep $UP_coord|awk '($3> '"$SS5"')')`
				if [[ ! -z "$mvp" ]]
				then
					echo $yy
					DN_coord=`echo $mvp|awk '{print $3}'`
					DN_count=0
					SK_count=`echo $mvp|awk '{print $4}'`
				elif [ -z "$mvp" ]
				then 
					DN_coord=0
					DN_count=0
					SK_count=0
				else
					echo "what?"
				fi
				samtools view -b /Users/rongzheng/Downloads/CryEx_2023/HTE1_chr21_sorted.bam "chr21:$SS3-$SS5"> GENE	
				Ex_count=`echo $(samtools view GENE|awk '(($4*1)>='"$SS3"' && ($4*1)<='"$SS5"')'|awk '($6=="100M")'|wc -l)`
				add_2=`echo "scale=10; ($UP_count+$Ex_count+$DN_count)/($SS5-$SS3+100)"|bc`
				tot_2=`echo "scale=10; $add_2+($SK_count/99)"|bc`
				PSI=`echo "scale=2; $add_2/$tot_2"|bc`
				printf "$chr\t$SS3\t$SS5\t$UP_coord\t$DN_coord\t$UP_count\t$DN_count\t$SK_count\t$Ex_count\t$PSI\t$sign\t$path\n" >> middle_first_minus.PSI
			done < UP
		elif [ -z "$max_SS3" ]
		then
			
			echo "noSS3only"
		else
			echo "what?"
		fi
				
		cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc3SS_HTE_count.txt > SS3_file
		while read ss;
		do	
			echo $ss
			SS3=`echo $ss|awk '{print $1}'`
			cat add_SS_count.txt|grep $SS3|awk '($2< '"$SS3"')' > UP
			while read ll;
			do
				echo $ll
				UP_coord=`echo $ll|awk '{print $2}'`
				UP_count=`echo $ll|awk '{print $4}'`
				cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc5SS_HTE_count.txt|awk '($1> '"$SS3"')' > SS5_file
				sud=`echo $(cat /Users/rongzheng/Downloads/CryEx_2023/Code/HTE1/sub/${path}final_junc5SS_HTE_count.txt|awk '($1> '"$SS3"')')`
				if [[ ! -z "$sud" ]]
				then
				
					while read mm;
					do
						echo $mm
						SS5=`echo $mm|awk '{print $1}'`
						cat add_SS_count.txt|grep $SS5|awk '($3> '"$SS5"')'> DN
						while read nn;
						do
							echo $nn
							DN_coord=`echo $nn|awk '{print $3}'`
							DN_count=`echo $nn|awk '{print $4}'`
							tmp=`echo $(cat add_SS_count.txt|awk '($2=='"$UP_coord"')'|awk '($3=='"$DN_coord"')')` 
					
							if [[ ! -z "$tmp" ]]
							then
								echo $tmp
								SK_count=`echo $tmp|awk '{print $4}'`
								echo "works"
							elif [ -z "$tmp" ]
							then
								SK_count=0
								echo "zero"
							else
								echo "what?"
							fi
							samtools view -b /Users/rongzheng/Downloads/CryEx_2023/HTE1_chr21_sorted.bam "chr19:$SS3-$SS5"> GENE	
							Ex_count=`echo $(samtools view GENE|awk '(($4*1)>='"$SS3"' && ($4*1)<='"$SS5"')'|awk '($6=="100M")'|wc -l)`
							add_2=`echo "scale=10; ($UP_count+$Ex_count+$DN_count)/($SS5-$SS3+100)"|bc`
							tot_2=`echo "scale=10; $add_2+($SK_count/99)"|bc`
							PSI=`echo "scale=2; $add_2/$tot_2"|bc`
						printf "$chr\t$SS3\t$SS5\t$UP_coord\t$DN_coord\t$UP_count\t$DN_count\t$SK_count\t$Ex_count\t$PSI\t$sign\t$path\n" >> middle_minus.PSI

						done < DN

					done < SS5_file
				elif [ -z "$sud" ]
				then
					echo "didbefore"
				else
					echo "what?"
				fi
				
			done < UP
			
		done < SS3_file
	else
		echo "what?"
	fi
done < trim_middle.coord
cat middle_minus.PSI middle_plus.PSI|awk '!seen[$0]++' > middle.PSI
cat middle_first_plus.PSI middle_fist_minus.PSI|awk '!seen[$0]++' > middle_first.PSI
cat middle_last_plus.PSI middle_last_plus.PSI|awk '!seen[$0]++' > middle_last.PSI
#rm middle_minus.PSI
#rm middle_plus.PSI


rm SK
rm UP
rm DN
rm chr
rm start
rm end
rm sign
rm lable





