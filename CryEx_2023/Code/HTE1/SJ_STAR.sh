#!/bin/bash

### SS_count.txt generates splice junction information for any expressed exon in each sample
samtools view ./CryEx_2023/HTE1_chr21_sorted.bam|awk '($6~/N/)'>spliced_reads.txt
cat spliced_reads.txt |awk '{print $4}' > position.txt
cat spliced_reads.txt |awk '{print $6}'|cut -f1 -d"M" > Match.txt
paste -d '\t' position.txt Match.txt > junc5SS.txt
cat spliced_reads.txt |awk '{print $6}'|cut -f1 -d"N"|sed 's/.*M//'>Skip.txt
paste -d '\t' junc5SS.txt Skip.txt > junc3SS.txt
cat junc3SS.txt |awk '{print $1+$2-1}'>junc5SS
cat junc3SS.txt |awk '{print $1+$2+$3-1}'>junc3SS
paste -d '\t' junc5SS junc3SS> junc.txt
cat junc.txt |awk '{count[$0]++} END {for (word in count) print word, count[word]}'|sort -k1,1n > junc_count.txt
cat junc_count.txt |awk '{print $1":"$2}' > SS
paste -d '\t' SS junc_count.txt > SS_count.txt
rm spliced_reads.txt
rm position.txt
rm Match.txt
rm junc5SS.txt
rm junc3SS.txt
rm Skip.txt
rm junc5SS
rm junc3SS
rm junc.txt
rm junc_count.txt
rm SS

### HTE1_SJ.out.tab is from STAR output containing splice junction information 
### SJ.tab contains chromosome coordinates and reads aligned to STAR predicted splice junction for chromosome 21
### comparing our generated any expressed exons' splice junction information with STAR output, we are able to retrieve an even more comprehensive cohort of splice junction information of all exons identified in a sample in add_SS_count.txt

cat ./CryEx_2023/HTE1_SJ.out.tab |awk '{print $1"\t"$2"\t"$3"\t"$7}'|awk '($1=="chr21")' > SJ.tab
cat SS_count.txt |awk '{print "chr21""\t"($2+1)"\t"$3"\t"$4}' > SS_count.tab
bedtools intersect -a SJ.tab -b SS_count.tab -F 1 -f 1 -v > tmp.SJ
cat tmp.SJ SS_count.tab |awk '{print "chr21""\t"($2-1)"\t"$3"\t"$4"\t"$2":"$3}' > add_SS_count.txt
rm SJ.tab
rm tmp.SJ
rm SS_count.tab
#rm SS_count.txt






