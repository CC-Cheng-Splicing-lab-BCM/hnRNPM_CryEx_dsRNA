##############################EMBO ssGSEA method##################################
library("htmltools")
#BiocManager::install("GSVA")
library("GSVA")

filename = "~/Downloads/ISG_manual.gmt.txt"
data = file(filename,open="r")
n = 1
pathway = list()
while(TRUE){
  line = readLines(data,n=1)
  if(length(line) == 0){
    break
  }
  pathway[n] = line
  n = n+1
}
close(data)
all_pathway = list()
names = c()
for(i in 1:length(pathway)){
  one_pathway = strsplit(as.character(pathway[i]),"\t")[[1]]
  names[i] = one_pathway[1]
  one_pathway = one_pathway[-1]
  all_pathway[[i]] = one_pathway
}
names(all_pathway) = names


# ################################################################################
# To take DLBC as an example in plotting ISG score of each sample in this cancer 
cancer_types = c('DLBC')
# define working directory
base_dir = './'

# counts table
counts_filename = sapply(cancer_types, function(x) paste('~/Downloads/', x, '_out_counts.txt', sep = ''))
# counts metadata
counts_meta = sapply(cancer_types, function(x) paste('~/Downloads/gdc_', x, '_counts.tsv', sep = ''))


###########################
i=1
DLBC_meta_ori<-read.table(counts_meta[[i]], sep = '\t', header = T)
colnames(DLBC_meta_ori)<-c('File.ID', 'File.Name', 'Data.Category', 'Data.Type', 'Project.ID', 
                          'Case.ID', 'Sample.ID', 'Sample.Type')
DLBC_meta<-DLBC_meta_ori[DLBC_meta_ori$Sample.Type=="Primary Tumor"|DLBC_meta_ori$Sample.Type=="Primary Blood Derived Cancer - Peripheral Blood",]
DLBC_meta$File.Name<-gsub("-", ".", substring(DLBC_meta$File.Name,1,36))
rownames(DLBC_meta)<-DLBC_meta$File.Name
DLBC_meta$Sample.ID<-substring(DLBC_meta$Sample.ID,1,15)
head(DLBC_meta)
dim(DLBC_meta)

# # ############################
# # # to extract Basal subtype of TCGA BRCA samples:
# metadata<-read.csv('~/Downloads/clinical_data_BRCA.txt', sep = '\t')
# # PAM50Call_RNAseq
# # sample_type
# unique(metadata$PAM50Call_RNAseq)
# LumA_ori<-metadata[metadata$PAM50Call_RNAseq=="Basal",c('sampleID','sample_type')]
# LumA<-LumA_ori[LumA_ori$sample_type=="Primary Tumor",]
# # get rid of a metastatic sample
# dim(LumA)
# 
# BRCA_meta_new<-merge(BRCA_meta, LumA, by.x="Sample.ID", by.y="sampleID")
# BRCA_meta_new<-BRCA_meta_new[!duplicated(BRCA_meta_new$Sample.ID),]
# dim(BRCA_meta_new)
# head(BRCA_meta_new)
################################################################################
################Clean up the counts table#######
DLBC_count<-read.table(counts_filename[[i]], sep = '\t', header = T)
DLBC_count$Ensembl_gene_id<-substring(DLBC_count$Ensembl_gene_id,1,15)
DLBC_count<-DLBC_count[!duplicated(DLBC_count$Ensembl_gene_id),]
rownames(DLBC_count)<-DLBC_count$Ensembl_gene_id
colnames(DLBC_count)<-gsub("X","", colnames(DLBC_count))

colnames(DLBC_count)<-c('Ensembl_gene_id', gsub(".htseq.counts","",colnames(DLBC_count[,-1])))
colnames(DLBC_count)<-c('Ensembl_gene_id', gsub("\\.","-",colnames(DLBC_count[,-1])))
dim(DLBC_count)
head(DLBC_count[1:3,1:3])
###########Extract Basal subtype table##########################
DLBC_meta_new$File.Name<-gsub('\\.', '-', DLBC_meta_new$File.Name)
DLBC_count_trim<-DLBC_count[,c(DLBC_meta_new$File.Name)]
head(DLBC_count_trim[,1:5])
# 60483
dim(DLBC_count_trim)
# ##################################################
#BiocManager::install('DESeq2')
library('DESeq2')
#BiocManager::install('edgeR')
library('edgeR')
#DLBC_count_trim<-risk_tot

input<-na.omit(DLBC_count_trim)

conditions<-na.omit(factor(colnames(DLBC_count_trim)))

colData<-na.omit(data.frame(sampleNames=colnames(input), conditions=conditions))
row.names(colData)<-colData$sampleNames
dds<-DESeqDataSetFromMatrix(countData = as.matrix(input),
                            colData = colData,
                            design= ~conditions)
y<-DGEList(counts(dds))
y$samples
# filter out low abundance genes
keep<-rowSums(cpm(counts(dds))>5/(min(y$samples$lib.size)/10^6))>(length(y$samples$lib.size)/2)
dds_kpt<-dds[keep,]
dim(dds_kpt)
##################lcpm TMM normalization################
ym<-DGEList(counts(dds_kpt))
ym<-calcNormFactors(ym, method="TMM")
lcpmtmm.rna<-cpm(ym, log=T, normalized.lib.sizes=T)
# hnRNPM: ENSG00000099783
# upper 33% quantile and lower 33% quantile of hnRNPM expression among TCGA DLBC samples
summary(as.numeric(lcpmtmm.rna[rownames(lcpmtmm.rna)=="ENSG00000099783",]))
low_q<-as.numeric(quantile(as.numeric(lcpmtmm.rna[rownames(lcpmtmm.rna)=="ENSG00000099783",]), probs = c(0.33, 0.67))[1])
high_q<-as.numeric(quantile(as.numeric(lcpmtmm.rna[rownames(lcpmtmm.rna)=="ENSG00000099783",]), probs = c(0.33, 0.67))[2])
lq<-c()
hq<-c()

for (i in 1:ncol(lcpmtmm.rna)){
  if (as.numeric(lcpmtmm.rna[rownames(lcpmtmm.rna)=="ENSG00000099783", i])< low_q){
    lq <- c(lq,i)
    
  } else if (as.numeric(lcpmtmm.rna[rownames(lcpmtmm.rna)=="ENSG00000099783", i])> high_q){
    hq <- c(hq, i)
  } else {
    print ("discard")
  }
}

length(lq)
length(hq)
lq_hq_df<-lcpmtmm.rna[,c(lq,hq)]
dim(lq_hq_df)
lq_hq_df[1:3,1:3]
lq_hq_label<-data.frame('File.Name'=colnames(lq_hq_df), 'label'=c(rep('low',length(lq)), rep('high', length(hq))))
head(lq_hq_label)
lq_hq_meta<-merge(lq_hq_label, DLBC_meta_new[,c(1,3)], by='File.Name')
lq_hq_meta<-lq_hq_meta[!duplicated(lq_hq_meta$File.Name),]
dim(lq_hq_meta)
head(lq_hq_meta)
################################################################################
### ISGs score calculation with ssGSEA
colnames(lq_hq_df)<-gsub("X","",colnames(lq_hq_df))
head(lq_hq_df[,1:5])
lcpmtmm<-data.frame(cbind('ensembl_gene_id'=rownames(lq_hq_df), lq_hq_df))
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),
  values=lcpmtmm$ensembl_gene_id,
  mart=mart)
head(genes)
dim(genes)
tmp_gene<-merge(lcpmtmm, genes, by='ensembl_gene_id')
tmp_gene<-tmp_gene[!duplicated(tmp_gene$hgnc_symbol),]
head(tmp_gene)
dim(tmp_gene)
gene<-tmp_gene[,-c(1,ncol(tmp_gene), (ncol(tmp_gene)-1))]
rownames(gene)<-tmp_gene$hgnc_symbol
colnames(gene)<-gsub("X","",colnames(gene))
head(gene[,1:4])
for (i in 1:ncol(gene)){
  gene[,i]<-as.numeric(gene[,i])
}
library(GSVA)
gsva_es <- gsva(as.matrix(gene), gset.idx.list=all_pathway, method="ssgsea",parallel.sz=4)
#########Using ISG scores to make boxplots##############################
ISGs_LUMA_input<-data.frame(cbind('ISGs_score'=as.numeric(gsva_es), 'sample.ID'=gsub("X","",colnames(gsva_es)),
                                  'cancers'=c(rep("DLBC_Basal", nrow(gsva_es))), 'hnM'=c(as.numeric(lcpmtmm[lcpmtmm$ensembl_gene_id=="ENSG00000099783",-1]))
))
ISGs_LUMA_input$ISGs_score<-as.numeric(ISGs_LUMA_input$ISGs_score)
ISGs_LUMA_input$hnM<-as.numeric(ISGs_LUMA_input$hnM)
summary(ISGs_LUMA_input)
ISGs_LUMA_input$sample.ID<-gsub('\\.','-',ISGs_LUMA_input$sample.ID)
head(ISGs_LUMA_input)
dim(ISGs_LUMA_input)
ISGs_meta<-merge(lq_hq_meta, ISGs_LUMA_input, by.x='File.Name', by.y='sample.ID')
ISGs_meta<-ISGs_meta[!duplicated(ISGs_meta$File.Name),]
head(ISGs_meta)
dim(ISGs_meta)
### Comparing the mean of the ISG scores of hnRNPM lowly expressed samples with that in hnRNPM highly expressed samples
mean(ISGs_meta[ISGs_meta$label=="low",]$ISGs_score)
mean(ISGs_meta[ISGs_meta$label=="high",]$ISGs_score)
wilcox.test(ISGs_meta[ISGs_meta$label=="low",]$ISGs_score, ISGs_meta[ISGs_meta$label=="high",]$ISGs_score)
################################################################################
### plot boxplots on ISGs scores in DLBC Basal 
library(ggplot2)
p<-ggplot(ISGs_meta, aes(x=label, y=ISGs_score))+
  geom_boxplot()
p + geom_jitter(shape=16, position=position_jitter(0.2))


