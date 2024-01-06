### HMLE##
# Gene Expression 
# load in epithelial samples gene expression data
HTE <- read.table('~/Downloads/HTE_HNRNPM_GeneCounts.txt', header = F, sep = '\t')
dim(HTE)
HTE$V1<-substr(HTE$V1, 1, 15)
head(HTE)
colnames(HTE)<-c('gene','len','HTE1','HTE2','HTE3','HTE4')
#HH_5<-read.table('~/Downloads/HNRNPM_RAW_Counts/HH-5ReadsPerGene.out.tab', sep = '\t')
#HH_5$V1<-substr(HH_5$V1,1,15)
#head(HH_5)
#colnames(HH_5)<-c('gene','HH_5','V3','V4')
#HH_6<-read.table('~/Downloads/HNRNPM_RAW_Counts/HH-6ReadsPerGene.out.tab', sep = '\t')
#HH_6$V1<-substr(HH_6$V1,1,15)
#head(HH_6)
#colnames(HH_6)<-c('gene','HH_6','V3','V4')
HH5<-read.table('~/Downloads/HH5ReadsPerGene.out.tab', sep = '\t')
HH5$V1<-substr(HH5$V1,1,15)
head(HH5)
colnames(HH5)<-c('gene','HH5','V3','V4')
HH6<-read.table('~/Downloads/HH6ReadsPerGene.out.tab', sep = '\t')
HH6$V1<-substr(HH6$V1,1,15)
head(HH6)
colnames(HH6)<-c('gene','HH6','V3','V4')
#####################################
#hnM_tot<-merge(merge(merge(merge(HTE, HH_5[,1:2], by='gene'), HH_6[,1:2], by='gene'), HH5[,1:2], by='gene'), HH6[,1:2], by='gene')
hnM_tot<-merge(merge(HTE,HH5[,1:2], by='gene'), HH6[,1:2], by='gene')
hnM_tot<-hnM_tot[!duplicated(hnM_tot$gene),]
rownames(hnM_tot)<-hnM_tot$gene
head(hnM_tot)
# 58609
dim(hnM_tot)
#### get rid of low abundance genes
keep <- rowSums(cpm(hnM_tot[,-c(1:2)])>0.5)>=3
hnM_keep <- hnM_tot[keep,]
head(hnM_keep)
# 15639
dim(hnM_keep)
kp_y<-DGEList(hnM_keep[,-c(1:2)])
kp_y$samples
kp_y<-calcNormFactors(kp_y, method='TMM')
lcpmtmm_hnM<-cpm(kp_y, log=F, normalized.lib.sizes=T)
dim(lcpmtmm_hnM)
head(lcpmtmm_hnM)
################################################
#sum(hnM_keep$HH_5)/1000000, sum(hnM_keep$HH_6)/1000000, 
sum_lib<-c(sum(hnM_keep$HTE1)/1000000, sum(hnM_keep$HTE2)/1000000, sum(hnM_keep$HTE3)/1000000, sum(hnM_keep$HTE4)/1000000,
           sum(hnM_keep$HH5)/1000000, sum(hnM_keep$HH6)/1000000)
hnM_norm<-hnM_keep
for (i in 3:ncol(hnM_norm)){
  hnM_norm[,i]<-hnM_keep[,i]/sum_lib[i-2]
}
hnM_norm_norm<-hnM_norm
for (k in 1:nrow(hnM_norm)){
  hnM_norm_norm[k,3:ncol(hnM_norm_norm)]<-as.numeric((hnM_norm[k,3:ncol(hnM_norm)]*1000))/hnM_norm$len[k]
}
#hnM_norm<-hnM_tot_mod[rowSums(hnM_tot_mod[,-c(1:2)]>5)>=1,]
dim(hnM_norm_norm)
head(hnM_norm_norm)
#  cat hallmark_interferon_alpha_response_gene_set.txt hallmark_interferon_gamma_response_gene_set.txt A.J.Scadden_nature_struct_bio.interferon_I_list|awk '!seen[$0]++' > ISG.txt
ISG<-read.table('~/Downloads/ISG.txt')
head(ISG)
hnM_ISG<-merge(symbol_HMLE, ISG, by.x='hgnc_symbol', by.y='V1')
hnM_ISG<-hnM_ISG[!duplicated(hnM_ISG$hgnc_symbol),]
head(hnM_ISG)
# 195
dim(hnM_ISG)
ISG_input<-hnM_ISG[,-c(1:3)]
rownames(ISG_input)<-hnM_ISG$hgnc_symbol
dim(ISG_input)
head(ISG_input)
################################################################################
################################################################################
################################################################################
for (i in 1:ncol(ISG_input)){
  ISG_input[,i]<-as.numeric(ISG_input[,i])
}
summary(ISG_input)
FC_input<-ISG_input
# 195; 9
dim(FC_input)
for (k in 1:nrow(FC_input)){
  FC_input[k,1]<-(ISG_input[k,1]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5,7)]))+0.01)
  FC_input[k,2]<-(ISG_input[k,2]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5,7)]))+0.01)
  FC_input[k,3]<-(ISG_input[k,3]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5,7)]))+0.01)
  FC_input[k,4]<-(ISG_input[k,4]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5,7)]))+0.01)
  FC_input[k,5]<-(ISG_input[k,5]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5,7)]))+0.01)
  FC_input[k,6]<-(ISG_input[k,6]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5,7)]))+0.01)
  FC_input[k,7]<-(ISG_input[k,7]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5,7)]))+0.01)
  FC_input[k,8]<-(ISG_input[k,8]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5,7)]))+0.01)
}
mean_KD<-c()
for (i in 1:nrow(FC_input)){
  mean_KD<-c(mean_KD, mean(as.numeric(FC_input[i,5:8])))
}
FC_input$mean_KD<-mean_KD
head(FC_input)
summary(FC_input)
new_rev_trim<-FC_input[rev(order(FC_input$mean_KD)),-ncol(FC_input)]
new_rev_trim<-na.omit(new_rev_trim)
new_rev_trim[new_rev_trim>=3]<-3
summary(new_rev_trim)
head(new_rev_trim)
rev_rev_trim<-new_rev_trim[rev(order(rowSums(new_rev_trim))),]
### ??? What is the method they use in Fig6C ??? 
library(pheatmap)
library(RColorBrewer)
#install.packages('viridis')
library(viridis)
genelist<-c("DDX58","BST2","IFI27","ADAR1","IFI6","IFIH1","IFIT1","IFIT2","IFIT3","IFITM1",
            "IRF7","ISG15","MX1","OAS3","OASL","STAT1")
labels<-rownames(rev_rev_trim)
labels[!labels %in% genelist] <- ""
#breaksList = seq(0, 100, by = 1)
# [1:50,]
pheatmap(mat=rev_rev_trim,
         #labels_row = labels[1:50],
         scale='none',
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         color=colorRampPalette(c("royalblue","white","red","red3"))(n=30),
         cluster_rows = FALSE,
         cluster_cols = TRUE)
################################################################################
# Try Quantile Normalization:
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

df<-hnM_norm_norm[,-c(1:2)]
rownames(df)<-hnM_norm_norm$gene
df<-lcpmtmm_hnM
head(df)
dim(df)
hnM_keep_df<-quantile_normalisation(df)
head(hnM_keep_df)
dim(hnM_keep_df)
summary(hnM_keep_df)
# convert Ensembl gene id to entrez ID
#library(biomaRt)
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),
  values=rownames(hnM_keep_df),
  mart=mart)
head(genes)
dim(genes)
symbol_HMLE<-merge(genes, data.frame(cbind('ensembl_gene_id'=rownames(hnM_keep_df), hnM_keep_df)), by="ensembl_gene_id")
symbol_HMLE<-symbol_HMLE[!duplicated(symbol_HMLE$ensembl_gene_id),]
head(symbol_HMLE)


#  cat hallmark_interferon_alpha_response_gene_set.txt hallmark_interferon_gamma_response_gene_set.txt A.J.Scadden_nature_struct_bio.interferon_I_list|awk '!seen[$0]++' > ISG.txt
ISG<-read.table('~/Downloads/ISG.txt')
head(ISG)
ISG_df<-merge(symbol_HMLE, ISG, by.x='hgnc_symbol', by.y='V1')
ISG_df<-ISG_df[!duplicated(ISG_df$hgnc_symbol),]
head(ISG_df)
# 195
dim(ISG_df)
ISG_input<-ISG_df
rownames(ISG_input)<-ISG_df$hgnc_symbol
dim(ISG_input)
head(ISG_input)
################################################################################
################################################################################
################################################################################
for (i in 4:ncol(ISG_input)){
  ISG_input[,i]<-as.numeric(ISG_input[,i])
}
summary(ISG_input)
FC_input<-ISG_input[,-c(1:3)]
# 195; 9
dim(FC_input)
head(FC_input)
for (k in 1:nrow(FC_input)){
  FC_input[k,1]<-(ISG_input[k,4]+0.01)/(mean(as.numeric(ISG_input[k,c(4,6,8)]))+0.01)
  FC_input[k,2]<-(ISG_input[k,5]+0.01)/(mean(as.numeric(ISG_input[k,c(4,6,8)]))+0.01)
  FC_input[k,3]<-(ISG_input[k,6]+0.01)/(mean(as.numeric(ISG_input[k,c(4,6,8)]))+0.01)
  FC_input[k,4]<-(ISG_input[k,7]+0.01)/(mean(as.numeric(ISG_input[k,c(4,6,8)]))+0.01)
  FC_input[k,5]<-(ISG_input[k,8]+0.01)/(mean(as.numeric(ISG_input[k,c(4,6,8)]))+0.01)
  FC_input[k,6]<-(ISG_input[k,9]+0.01)/(mean(as.numeric(ISG_input[k,c(4,6,8)]))+0.01)
}

mean_KD<-c()
for (i in 1:nrow(FC_input)){
  mean_KD<-c(mean_KD, mean(as.numeric(FC_input[i,c(2,4,6)])))
}
FC_input$mean_KD<-mean_KD
head(FC_input)
summary(FC_input)
new_rev_trim<-FC_input[rev(order(FC_input$mean_KD)),-ncol(FC_input)]
new_rev_trim<-na.omit(new_rev_trim)
summary(new_rev_trim)
new_rev_trim[new_rev_trim>=3]<-3
summary(new_rev_trim)
head(new_rev_trim)
rev_rev_trim<-new_rev_trim[rev(order(rowSums(new_rev_trim))),]
### ??? What is the method they use in Fig6C ??? 
library(pheatmap)
library(RColorBrewer)
#install.packages('viridis')
library(viridis)
genelist<-c("DDX58","BST2","IFI27","ADAR1","IFI6","IFIH1","IFIT1","IFIT2","IFIT3","IFITM1",
            "IRF7","ISG15","MX1","OAS3","OASL","STAT1")
labels<-rownames(rev_rev_trim)
labels[!labels %in% genelist] <- ""
#breaksList = seq(0, 100, by = 1)
# [1:50,]
pheatmap(mat=rev_rev_trim,
         #labels_row = labels[1:50],
         scale='none',
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         color=colorRampPalette(c("royalblue","white","red","red3"))(n=30),
         cluster_rows = FALSE,
         cluster_cols = TRUE)
################################################################################



