##############################EMBO##################################
library("htmltools")
#BiocManager::install("GSVA")
library("GSVA")
# output is empty
#myoutf <- "~/Downloads/ACC_score.xls"
for(i in 1 : length(all_pathway))
{
  cat("\r",i)
  com <- intersect(as.vector(unlist(all_pathway[[i]])), gene)
  all_pathway[[i]] <- com
}


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
cancer_types = c('UVM')
# define working directory
base_dir = '~/Downloads/'

# counts table
counts_filename = sapply(cancer_types, function(x) paste(base_dir, 'GDC_', x,'_counts/',x,'_out_counts.txt',sep = ''))                                 
# counts metadata
counts_meta = sapply(cancer_types, function(x) paste('~/Downloads/gdc_', x, '_counts.tsv', sep = ''))


###########################
i=1
ACC_meta_ori<-read.table(counts_meta[[i]], sep = '\t', header = T)
colnames(ACC_meta_ori)<-c('File.ID', 'File.Name', 'Data.Category', 'Data.Type', 'Project.ID', 
                          'Case.ID', 'Sample.ID', 'Sample.Type')
ACC_meta<-ACC_meta_ori[ACC_meta_ori$Sample.Type=="Primary Tumor"|ACC_meta_ori$Sample.Type=="Primary Blood Derived Cancer - Peripheral Blood",]
ACC_meta$File.Name<-gsub("-", ".", substring(ACC_meta$File.Name,1,36))
rownames(ACC_meta)<-ACC_meta$File.Name
ACC_meta$Sample.ID<-substring(ACC_meta$Sample.ID,1,15)
head(ACC_meta)
dim(ACC_meta)

# # ############################
# # # to extract normal samples:
# metadata<-read.csv('~/Downloads/clinical_data_BRCA.txt', sep = '\t')
# # PAM50Call_RNAseq
# # sample_type
# unique(metadata$PAM50Call_RNAseq)
# LumA_ori<-metadata[metadata$PAM50Call_RNAseq=="Basal",c('sampleID','sample_type')]
# LumA<-LumA_ori[LumA_ori$sample_type=="Primary Tumor",]
# # get rid of a metastatic sample
# dim(LumA)
# 
# ACC_meta_new<-merge(ACC_meta, LumA, by.x="Sample.ID", by.y="sampleID")
# ACC_meta_new<-ACC_meta_new[!duplicated(ACC_meta_new$Sample.ID),]
# dim(ACC_meta_new)
# head(ACC_meta_new)
################################################################################
######lcpmtmm normalization#####
############################
ACC_count<-read.table(counts_filename[[i]], sep = '\t', header = T)
ACC_count$Ensembl_gene_id<-substring(ACC_count$Ensembl_gene_id,1,15)
ACC_count<-ACC_count[!duplicated(ACC_count$Ensembl_gene_id),]
rownames(ACC_count)<-ACC_count$Ensembl_gene_id
colnames(ACC_count)<-gsub("X","", colnames(ACC_count))
#### get rid of low abundance genes
colnames(ACC_count)<-c('Ensembl_gene_id', gsub(".htseq.counts","",colnames(ACC_count[,-1])))
# # CPTAC
colnames(ACC_count)<-c('Ensembl_gene_id', gsub("\\.","-",colnames(ACC_count[,-1])))
dim(ACC_count)
head(ACC_count[1:3,1:3])
###########=====================================================================
ACC_meta$File.Name<-gsub('\\.', '-', ACC_meta$File.Name)
ACC_count_trim<-ACC_count[,c(ACC_meta$File.Name)]
# ACC_meta_new$File.Name<-gsub('\\.', '-', ACC_meta_new$File.Name)
# ACC_count_trim<-ACC_count[,c(ACC_meta_new$File.Name)]
head(ACC_count_trim[,1:5])
# 60483
dim(ACC_count_trim)
# ##################################################
#BiocManager::install('DESeq2')
#library('DESeq2')
#BiocManager::install('edgeR')
#library('edgeR')
#ACC_count_trim<-risk_tot

input<-na.omit(ACC_count_trim)

conditions<-na.omit(factor(colnames(ACC_count_trim)))

colData<-na.omit(data.frame(sampleNames=colnames(input), conditions=conditions))
row.names(colData)<-colData$sampleNames
dds<-DESeqDataSetFromMatrix(countData = as.matrix(input),
                            colData = colData,
                            design= ~conditions)
y<-DGEList(counts(dds))
y$samples
keep<-rowSums(cpm(counts(dds))>5/(min(y$samples$lib.size)/10^6))>(length(y$samples$lib.size)/2)
dds_kpt<-dds[keep,]
# 24592
dim(dds_kpt)
# #################
ym<-DGEList(counts(dds_kpt))
ym<-calcNormFactors(ym, method="TMM")
lcpmtmm.rna<-cpm(ym, log=T, normalized.lib.sizes=T)
#write.table(lcpmtmm.rna, '~/Downloads/BRCA_Basal_lcpmtmm.rna.txt', sep = '\t', row.names = T, quote = F)
# hnRNPM: ENSG00000099783
summary(as.numeric(lcpmtmm.rna[rownames(lcpmtmm.rna)=="ENSG00000099783",]))
low_q<-as.numeric(quantile(as.numeric(lcpmtmm.rna[rownames(lcpmtmm.rna)=="ENSG00000099783",]), probs = c(0.25, 0.75))[1])
high_q<-as.numeric(quantile(as.numeric(lcpmtmm.rna[rownames(lcpmtmm.rna)=="ENSG00000099783",]), probs = c(0.25, 0.75))[2])
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
# 35
length(lq)
# 35
length(hq)
lq_hq_df<-lcpmtmm.rna[,c(lq,hq)]
dim(lq_hq_df)
lq_hq_df[1:3,1:3]
lq_hq_label<-data.frame('File.Name'=colnames(lq_hq_df), 'label'=c(rep('low',length(lq)), rep('high', length(hq))))
head(lq_hq_label)
#lq_hq_meta<-merge(lq_hq_label, ACC_meta_new[,c(1,3)], by='File.Name')
lq_hq_meta<-merge(lq_hq_label, ACC_meta[,c(2,7)], by='File.Name')
lq_hq_meta<-lq_hq_meta[!duplicated(lq_hq_meta$File.Name),]
dim(lq_hq_meta)
head(lq_hq_meta)
################################################################################
### ISGs 
#lcpmtmm.rna<-read.table('~/Downloads/BRCA_Basal_lcpmtmm.rna.txt', sep = '\t', header = T)
colnames(lq_hq_df)<-gsub("X","",colnames(lq_hq_df))
head(lq_hq_df[,1:5])
lcpmtmm<-data.frame(cbind('ensembl_gene_id'=rownames(lq_hq_df), lq_hq_df))
#head(lcpmtmm)
#library(biomaRt)
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
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
##############################################################
#library(GSVA)
gsva_es <- gsva(as.matrix(gene), gset.idx.list=all_pathway, method="ssgsea",parallel.sz=4)
ISGs_LUMA_input<-data.frame(cbind('ISGs_score'=as.numeric(gsva_es), 'sample.ID'=gsub("X","",colnames(gsva_es)),
                                 'cancers'=c(rep("UVM", nrow(gsva_es))), 'hnM'=c(as.numeric(lcpmtmm[lcpmtmm$ensembl_gene_id=="ENSG00000099783",-1]))
))
ISGs_LUMA_input$ISGs_score<-as.numeric(ISGs_LUMA_input$ISGs_score)
ISGs_LUMA_input$hnM<-as.numeric(ISGs_LUMA_input$hnM)
summary(ISGs_LUMA_input)
ISGs_LUMA_input$sample.ID<-gsub('\\.','-',ISGs_LUMA_input$sample.ID)
head(ISGs_LUMA_input)
dim(ISGs_LUMA_input)
head(lq_hq_meta)
ISGs_meta<-merge(lq_hq_meta, ISGs_LUMA_input, by.x='File.Name', by.y='sample.ID')
ISGs_meta<-ISGs_meta[!duplicated(ISGs_meta$File.Name),]
head(ISGs_meta)
dim(ISGs_meta)

mean(ISGs_meta[ISGs_meta$label=="low",]$ISGs_score)
mean(ISGs_meta[ISGs_meta$label=="high",]$ISGs_score)
wilcox.test(ISGs_meta[ISGs_meta$label=="low",]$ISGs_score, ISGs_meta[ISGs_meta$label=="high",]$ISGs_score)

write.table(ISGs_meta, '~/Downloads/UVM_ISGs_score.txt', sep = '\t', row.names = F, quote = F)
################################################################################
### plot boxplots on ISGs across cancers
#library(ggplot2)
p<-ggplot(ISGs_meta, aes(x=label, y=ISGs_score))+
  geom_boxplot()
p + geom_jitter(shape=16, position=position_jitter(0.2))
# BRCA_Basal_ISGs_boxplot_70.pdf
##############################################################################################################################################
# Figure6B:
ISGs_cancer<-read.table('~/Downloads/BRCA_Basal_ISGs_score.txt', sep = '\t', header = T)
head(ISGs_cancer)
dim(ISGs_cancer)


library(ggplot2)
ggplot(ISGs_cancer, aes(x=hnM, y=EMT_score))+
  geom_smooth(method = lm, color='#00CCFF', se=TRUE)+
  geom_point(color='magenta', size=4, alpha=1)
#xlim(-0.5,1)+ylim(-0.5,1)
fit <- lm(EMT_score ~ hnM, ISGs_cancer)
summary(fit)
cor.test(ISGs_cancer$hnM, ISGs_cancer$EMT_score)

ggplot(ISG_df, aes(x=hnM_lb, y=EMT_score)) +
  geom_boxplot(notch = TRUE)
wilcox.test(ISG_df[ISG_df$hnM_lb=="low",]$EMT_score, ISG_df[ISG_df$hnM_lb=="high",]$EMT_score)

