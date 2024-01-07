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
## Output format
### 
## IGV and PCR validation
<img width="530" alt="Screen Shot 2024-01-07 at 12 02 09 AM" src="https://github.com/CC-Cheng-Splicing-lab-BCM/hnRNPM_CryEx_dsRNA/assets/45469780/37273376-d2dd-4d4a-8123-925148928a61">
<img width="412" alt="Screen Shot 2024-01-07 at 12 16 36 AM" src="https://github.com/CC-Cheng-Splicing-lab-BCM/hnRNPM_CryEx_dsRNA/assets/45469780/beb66a3e-6353-4473-bce5-39a716d1885c">

Part II: visualization analyses

### Thank you
Thank you for using CryEx pipeline and the visualization analyses!

### A python version of the CryEx pipeline that is more user-friendly is underdevelopment, please stay tuned! If you have any questions, please do not hesitate to contact me through GitHub or at my email: Rong.Zheng@bcm.edu
=======
This repository contains the novel CryEx pipeline that processes cryptic splicing analysis and visualization analyses on the hnRNPM-repressed CryEx study
>>>>>>> parent of 0f055419 (Update README.md)
