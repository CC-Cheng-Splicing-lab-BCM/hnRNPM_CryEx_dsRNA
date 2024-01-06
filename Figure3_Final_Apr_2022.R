################################################################################
### Stacked barplot: Fig3C ###
df<-t(data.frame(c(338/344, 6/344)))
colnames(df)<-c('deep','prox')
rownames(df)<-c('CryEx')
spineplot(df)

counts<-c(6, 338)
barplot(counts)
################################################################################
################################################################################
# length distribution of CryEx: (plus AS)
############AS###########
# cat annotated_Ctrl_hnM_KD.txt|awk '{print $4"\t"$6"\t"$7"\t"$2"\t"$3"\t"$5}'|awk '!seen[$0]++'> annotated_Ctrl_hnM_KD_exon_coord.bed
AS_hnM<-read.table('~/Downloads/annotated_Ctrl_hnM_KD_exon_coord.bed', sep = '\t')
head(AS_hnM)
dim(AS_hnM)
summary(AS_hnM)

AS_dist<-c()
for (z in 1:nrow(AS_hnM)){
  AS_dist<-c(AS_dist, (AS_hnM$V3[z] - AS_hnM$V2[z] + 1))
}
# 159.9853
mean(AS_dist)
################################################################################
################################################################################
last_dist<-c()
first_dist<-c()
middle_dist<-c()

CryEx_M<-read.table('~/Downloads/trim_check_CryEx_M.txt', sep = '\t')
head(CryEx_M)
dim(CryEx_M)

CryEx_L<-read.table('~/Downloads/trim_check_CryEx_L.txt', sep = '\t')
head(CryEx_L)
dim(CryEx_L)

CryEx_F<-read.table('~/Downloads/trim_check_CryEx_F.txt', sep = '\t')
head(CryEx_F)
dim(CryEx_F)

#gg<-data.frame('chr'=CryEx_L$V1, 'start'=CryEx_L$V2, 'end'=CryEx_L$V3, 'len'=(CryEx_L$V3-CryEx_L$V2+1))
#write.table(gg, '~/Downloads/gg.txt', sep = '\t')

for (i in 1:nrow(CryEx_M)){
  middle_dist<-c(middle_dist, (CryEx_M$V3[i]-CryEx_M$V2[i]+1))
}
for (k in 1:nrow(CryEx_L)){
  last_dist<-c(last_dist, (CryEx_L$V3[k]-CryEx_L$V2[k]+1))
}
for (m in 1:nrow(CryEx_F)){
  first_dist<-c(first_dist, (CryEx_F$V3[m]-CryEx_F$V2[m]+1))
}
# 129.5235
mean(middle_dist)
# 829.9104
mean(last_dist)
# 128.5208
mean(first_dist)
# 324.9479
mean(c(middle_dist, last_dist, first_dist))
################################################################################
input_len<-data.frame('types'=c(rep('middle', length(middle_dist)), rep('last', length(last_dist)),
                                rep('first', length(first_dist))),
                      'length'=log2(c(middle_dist, last_dist, first_dist)))
head(input_len)
dim(input_len)
summary(input_len)
### length distribution ###
#install.packages("ggstatsplot")
#library(ggstatsplot)
#install.packages("palmerpenguins")
#library(palmerpenguins)
#install.packages("tidyverse")
#library(tidyverse)
#install.packages("PMCMRplus")
#library(PMCMRplus)

#data("penguins", package = "palmerpenguins")

#penguins<-drop_na(penguins)
#penguins<-input_len
plt<-ggbetweenstats(
  data=input_len,
  x=types,
  y=length,
  p.adjust.method = "fdr"
)

plt<-plt+
  # add labels and title
  labs(
    x="Cryptic Exon Type",
    y="Exon Length",
    title = "Distribution of exon length across cryptic exon types"
  )
#+
# Customizations
#theme(
# This is the new default font in the plot
#family = "Arial",
#text = element_text(size = 8, color = "black"),
# plot.title = element_text(
#family = "Helvetica",
# size = 20,
# face = "bold",
# color = "#2a475e"
# ),
# Statistical annotations below the main title
# plot.subtitle = element_text(
#family = "Arial",
#  size = 15,
# face = "bold",
# color = "#1b2838"
#  ),
#  plot.title.position = "plot", # slightly different from default
# axis.text = element_text(size = 10, color = "black"),
# axis.title = element_text(size = 12)

# )

# 1. Remove axis ticks
# 2. Change default color of the axis lines with a lighter one
# 3. Remove most reference lines, only keep the major horizontal ones
#    This reduces clutter, while keeping the reference for the variable
#    being compared.
# 4. Set the panel and the background fill to the same light color.
plt<-plt+
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )
plt
ggsave(
  filename = here::here("Downloads", "tot_number_intron.pdf"),
  plot = plt,
  #width = 8,
  #height = 8,
  #device = "pdf"
)



################################################################################
###CD44 KO
LM2_fc<-read.table('~/Downloads/CCheng_LM2_fc.txt', sep = '\t', header = F)
LM2_fc$V1<-substring(LM2_fc$V1, 1, 15)
head(LM2_fc)
dim(LM2_fc)
LM2_fc<-LM2_fc[!duplicated(LM2_fc$V1),]
LM2_fc[LM2_fc$V1=="ENSG00000026508",]
colnames(LM2_fc)<-c('gene', 'Ctrl_rep1','KO_rep1','Ctrl_rep2','KO_rep2')

#BiocManager::install('DESeq2')
library('DESeq2')
#BiocManager::install('edgeR')
library('edgeR')


input_LM2<-LM2_fc[,-1]
rownames(input_LM2)<-LM2_fc$V1
colnames(input_LM2)<-c('Ctrl_rep1','KO_rep1','Ctrl_rep2','KO_rep2')

conditions<-factor(c('Ctrl','KO','Ctrl','KO'))

colData<-data.frame(sampleNames=colnames(input_LM2), conditions=conditions)

row.names(colData)<-colData$sampleNames
dds<-DESeqDataSetFromMatrix(countData = as.matrix(input_LM2),
                            colData = colData,
                            design= ~conditions)
y<-DGEList(counts(dds))
y$samples
keep<-rowSums(cpm(counts(dds))>5/(min(y$samples$lib.size)/10^6))>=2
dds_kpt<-dds[keep,]
# 17647
dim(dds_kpt)
#####################
dds_fin<-DESeq(dds_kpt)
res<-results(dds_fin)
head(res)
head(counts(dds_kpt))
res<-na.omit(res)
dim(res)
dds_count<-merge(LM2_fc, data.frame(cbind('gene'=rownames(res), res)), by='gene')
dds_count<-dds_count[!duplicated(dds_count$gene),]
head(dds_count)
dim(dds_count)
dim(dds_count[abs(dds_count$log2FoldChange)>1 & dds_count$padj<0.05,])
write.table(dds_count, '~/Downloads/dds_CD44_KO_Ctrl_LM2_Apr_7_2022.txt', sep = '\t', row.names = F)


### Add gene symbol and refseq
# library(biomaRt)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),
  values=dds_count$gene,
  mart=mart)
head(genes)
dim(genes)

dds_count_symbol<-merge(genes[,c(1,3)], dds_count, by.x='ensembl_gene_id', by.y='gene')
dds_count_symbol<-dds_count_symbol[!duplicated(dds_count_symbol$ensembl_gene_id),]
# 17624
dim(dds_count_symbol)
head(dds_count_symbol)

### Add refseq
refseq_summary<-read.table('~/Downloads/hg19_refseq_summaries.txt', sep = '\t', header = F)
head(refseq_summary[,1:2])
trim_refseq<-refseq_summary[!duplicated(refseq_summary$V2),]


refseq<-read.table("~/Downloads/refseq_nonrepeats.txt", header = T, na.strings = c("","NA"))
head(refseq[,1])
dim(refseq)


# 27939
colnames(refseq)<-c('hgnc_symbol','refseq_func')
dds_count_symbol_refseq<-merge(dds_count_symbol, refseq, by='hgnc_symbol', all=T)
dds_count_symbol_refseq<-dds_count_symbol_refseq[!duplicated(dds_count_symbol_refseq$ensembl_gene_id),]
dim(dds_count_symbol_refseq)
write.table(dds_count_symbol_refseq, 'LM2_refseq_CD44_KO_04072022.txt', sep = '\t', row.names = F)

################################################################################
################################################################################
################################################################################
# Adapted from Ule Cell 2018 paper:
# For the purpose of extracting repeat element chromosome coordinates from hg19 (UCSC)
setwd('~/Downloads/ULE_2020/')
system("pwd")
system("wget http://repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz")
system("gzip -d hg19.fa.out.gz")

library("data.table")
repeats.dt <- fread("hg19.fa.out", fill=TRUE, skip = 3L, header=FALSE)

#badly formatted table
repeats.dt <- repeats.dt [ 3:nrow(repeats.dt), ]
colnames(repeats.dt) <- c("swScore", "milliDiv", "milliDel", "milliIns", "seqnames", "start", "end", "geneLeft", "strand", "repName", "repFamily", "repStart", "repEnd", "repLeft", "bin")

repeats.dt [, repClass := tstrsplit(repFamily, "/", keep =1) ] 
# 5467455 repeats
repeats.dt [ strand == "C", strand := "-" ]
head(repeats.dt)
dim(repeats.dt)
write.table(repeats.dt, '~/Downloads/Repeats.bed', row.names = FALSE, col.names = FALSE)
#################################################################
#repeatsofInterest.dt <- repeats.dt [ repClass %in% c("DNA") , ]
repeatsofInterest.dt <- repeats.dt
dim(repeatsofInterest.dt)
repeatsofInterest.dt <- repeatsofInterest.dt [ !grepl("tRNA", repFamily), ]
# 1569213
dim(repeatsofInterest.dt)
#######################Figure 1A##################################
head(repeatsofInterest.dt)
classic_repeatsofInterest.dt<-rbind(repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr1",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr2",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr3",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr4",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr5",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr6",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr7",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr8",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr9",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr10",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr11",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr12",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr13",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr14",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr15",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr16",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr17",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr18",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr19",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr20",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr21",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chr22",],
                                    repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chrX",], repeatsofInterest.dt[repeatsofInterest.dt$seqnames=="chrY",])
head(classic_repeatsofInterest.dt)
# LINE: 1798202
# SINE: 1774227
# LTR: 732108
# DNA: 519257
# repeats: 5388919
dim(classic_repeatsofInterest.dt)
write.table(classic_repeatsofInterest.dt, '~/Downloads/LINEs.bed', row.names = FALSE, col.names = FALSE)
write.table(classic_repeatsofInterest.dt, '~/Downloads/SINEs.bed', row.names = FALSE, col.names = FALSE)
write.table(classic_repeatsofInterest.dt, '~/Downloads/LTR.bed', row.names = FALSE, col.names = FALSE)
write.table(classic_repeatsofInterest.dt, '~/Downloads/DNA.bed', row.names = FALSE, col.names = FALSE)
write.table(classic_repeatsofInterest.dt, '~/Downloads/Repeats_all.bed', row.names = FALSE, col.names = FALSE)
#####################################################################################################################
# cat LINEs.bed|awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}'|tr -d '"'> LINE_coord.bed
# cat SINEs.bed|awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}'|tr -d '"'> SINE_coord.bed
# cat LTR.bed|awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}'|tr -d '"'> LTR_coord.bed
# cat DNA.bed|awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}'|tr -d '"'> DNA_coord.bed
# cat Repeats_all.bed|awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}'|tr -d '"'> Repeats_coord.bed 
### Get deep intronic repeats:
# bedtools intersect -a gencode.v19_distal.intron -b LINE_coord.bed|awk '!seen[$0]++'> LINE_deep.bed
# bedtools intersect -a gencode.v19_distal.intron -b SINE_coord.bed|awk '!seen[$0]++'> SINE_deep.bed
# bedtools intersect -a gencode.v19_distal.intron -b LTR_coord.bed|awk '!seen[$0]++'> LTR_deep.bed 
# bedtools intersect -a gencode.v19_distal.intron -b DNA_coord.bed|awk '!seen[$0]++'> DNA_deep.bed
### Calculate GC content in deep intronic repeats:
















