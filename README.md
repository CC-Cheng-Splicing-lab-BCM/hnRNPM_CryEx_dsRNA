# CryEx_hnRNPM
<<<<<<< HEAD
This repository contains the novel bioinformatics pipeline named "CryEx" that identifies any expressed exon from RNA-seq and some visualization analyses on the hnRNPM-repressed CryEx study. 

Part I: CryEx

## Toy example dataset 
Example of RNA-seq data on chromosome 21 is available at GitHub HTE1_chr21_sorted.bam, HTE1_chr21_sorted.bam.bai

## Preparing the reference genome sequences and gene annotation files
### Human reference genome primary assembly (GRCh37)
### Human transcriptome GENCODE version 24 backmap 37 comprehensive assembly
## Run CryEx pipeline step by step
### Step 1: this step generates a "sub" folder that contains individual subfolder for each expressed exon identified
bash CryEx_Boundaries.sh
### Step 2: this step generates accurate 5' and 3' splice sites for any expressed exon
bash CryEx_Coor_Exact.sh
### Step 3: this step generates number of splice junctions in any expressed exon
bash CryEx_len.sh
### Step 3: this step generates splice junction information for any expressed exon (requires STAR output file SJ.out.tab)
bash SJ_STAR.sh
### Step 4: this step calculates PSI values for each expressed exon
bash SIMPLE_JUNC.sh
### Step 5: this step outputs PSI values for any identified intron including retained intron (requires intron chromosome coordinates file: SS_count.txt)
bash IR.sh
## Output format: 
## (chromosome/exon_start/exon_end/upstream_5'ss/downstream_3'ss/exon_reads/3'ss_spliced_read/5'ss_spliced_reads/skipped_reads/PSI/strandness)
### middle.PSI: contains any expressed exon that has both 5' and 3' splice sites
### first.PSI: contains any expressed exon that has 5' splice site only
### last.PSI: contains any expressed exon that has 3' splice site only
## (coordinate_ID/left_coordinate/right_coordinate/skipped_reads/intron_reads/intron_length/PSI)
### (chr21_)intron.PSI: contains any identifiable intron
## IGV and PCR validation
[IGV_RT_PCR.pdf](https://github.com/CC-Cheng-Splicing-lab-BCM/hnRNPM_CryEx_dsRNA/files/13852485/IGV_RT_PCR.pdf)

Part II: visualization analyses

### Thank you
Thank you for using CryEx pipeline and the visualization analyses!

### A python version of the CryEx pipeline that is more user-friendly is underdevelopment, please stay tuned! If you have any questions, please do not hesitate to contact me through GitHub or at my email: Rong.Zheng@bcm.edu
=======
This repository contains the novel CryEx pipeline that processes cryptic splicing analysis and visualization analyses on the hnRNPM-repressed CryEx study
>>>>>>> parent of 0f055419 (Update README.md)
