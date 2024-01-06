#!/bin/bash

### This script generate individual file that can obtain 5'SS and 3'SS junction reads and junction chromosome coordinates for each expressed exon in the sample
### This script is nested with CryEx_Boundaries.sh
coord=AXX
up_coord=UXX
down_coord=OXX
datadir=DXX

mkdir -p $datadir
cd $datadir
echo $coord

### Get Bam reads for each expressed exon region:
samtools view -b ./CryEx_2023/HTE1_chr21_sorted.bam ${coord}>HTE
### focus on reads that are spliced in this region:
samtools view HTE|awk '($6~/N/)' > HTE_spliced_reads.txt

### get the chromosome coordiates of the 5'SS junction and its reads number aligned there
cat HTE_spliced_reads.txt|awk '{print $4}'>position.txt
cat HTE_spliced_reads.txt|awk '{print $6}'|cut -f1 -d"M">Match.txt
paste -d '\t' position.txt Match.txt> junc5SS_HTE.txt
cat junc5SS_HTE.txt|awk '{print $1+$2-1}'|awk '{count[$1]++} END {for (word in count) print word, count[word]}'|sort -k2,2n> junc5SS_HTE_count.txt
cat junc5SS_HTE_count.txt|awk '($1 >= ('${up_coord}'-2) && $1 <= ('${down_coord}'+2))' > final_junc5SS_HTE_count.txt

### get the chromosome coordiates of the 3'SS junction and its reads number aligned there
cat HTE_spliced_reads.txt|awk '{print $6}'|cut -f1 -d"N"|sed 's/.*M//'>Skip.txt
paste -d '\t' junc5SS_HTE.txt Skip.txt > junc3SS_HTE.txt
cat junc3SS_HTE.txt|awk '{print $1+$2+$3-1}'|awk '{count[$1]++} END {for (word in count) print word, count[word]}'|sort -k2,2n>junc3SS_HTE_count.txt
cat junc3SS_HTE_count.txt|awk '($1 >= ('${up_coord}'-2) && $1 <= ('${down_coord}'+2))' > final_junc3SS_HTE_count.txt


echo 'junc5SS'

tail final_junc5SS_HTE_count.txt

echo 'junc3SS'

tail final_junc3SS_HTE_count.txt


exit 0

