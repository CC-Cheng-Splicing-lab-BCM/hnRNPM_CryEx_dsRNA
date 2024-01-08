### HMLE##
# Gene Expression 
# load in epithelial samples gene expression data
HTE <- read.table('~/Downloads/HTE_HNRNPM_GeneCounts.txt', header = F, sep = '\t')
dim(HTE)
HTE$V1<-substr(HTE$V1, 1, 15)
head(HTE)
colnames(HTE)<-c('gene','len','HTE1','HTE2','HTE3','HTE4')
HH5<-read.table('~/Downloads/HH5ReadsPerGene.out.tab', sep = '\t')
HH5$V1<-substr(HH5$V1,1,15)
head(HH5)
colnames(HH5)<-c('gene','HH5','V3','V4')
HH6<-read.table('~/Downloads/HH6ReadsPerGene.out.tab', sep = '\t')
HH6$V1<-substr(HH6$V1,1,15)
head(HH6)
colnames(HH6)<-c('gene','HH6','V3','V4')
#####################################
hnM_tot<-merge(merge(HTE,HH5[,1:2], by='gene'), HH6[,1:2], by='gene')
hnM_tot<-hnM_tot[!duplicated(hnM_tot$gene),]
rownames(hnM_tot)<-hnM_tot$gene
head(hnM_tot)
# 58609
dim(hnM_tot)

# batch effect correction
batch = c(1,1,1,1,2,2)
pheno<-data.frame('sample'=c(1,2,3,4,5,6), 'batch'=batch, 'cancer'=c('Ctrl','KD','Ctrl','KD','Ctrl','KD'))
rownames(pheno)<-c('HTE1','HTE2','HTE3','HTE4','HH5','HH6')
conditions = pheno$cancer
library_methods = pheno$batch
groups = sapply(as.character(conditions), switch, "Ctrl"=1, "KD"=2, USE.NAMES = F)
batches = as.numeric(pheno$batch)
library('sva')
#BiocManager::install('sva')
corrected_data = ComBat_seq(counts = as.matrix(hnM_tot[,-c(1:2)]),
                            batch = batches, group = groups)
head(corrected_data)
rownames(corrected_data)<-hnM_tot$gene

corrected_data<-data.frame(corrected_data)
head(corrected_data)
# fpkm normalization
sum_lib<-c(sum(corrected_data$HTE1), sum(corrected_data$HTE2), sum(corrected_data$HTE3), sum(corrected_data$HTE4), sum(corrected_data$HH5), sum(corrected_data$HH6))
#(5*1000000/(min(sum_lib)))
keep<-rowSums(cpm(corrected_data)>0.5)>=3
cor_df<-corrected_data[keep,]
head(cor_df)
dim(cor_df)
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),
  values=rownames(cor_df),
  mart=mart)
head(genes)
dim(genes)
symbol_HMLE<-merge(genes, data.frame("gene"=rownames(cor_df), cor_df), by.x="ensembl_gene_id", by.y="gene")
symbol_HMLE<-symbol_HMLE[!duplicated(symbol_HMLE$ensembl_gene_id),]
head(symbol_HMLE)
dim(symbol_HMLE)
symbol_HMLE<-symbol_HMLE[!duplicated(symbol_HMLE$hgnc_symbol),]
sum_HMLE<-c(sum(symbol_HMLE$HTE1)/1000000, sum(symbol_HMLE$HTE2)/1000000, sum(symbol_HMLE$HTE3)/1000000,
            sum(symbol_HMLE$HTE4)/1000000, sum(symbol_HMLE$HH5)/1000000, sum(symbol_HMLE$HH6)/1000000)
# divide the read counts by the "per million" scaling factor. 
HMLE_df<-symbol_HMLE[,-c(1:3)]
for (i in 1:ncol(HMLE_df)){
  HMLE_df[,i]<-symbol_HMLE[,(i+3)]/sum_HMLE[i]
}
rownames(HMLE_df)<-symbol_HMLE$hgnc_symbol
head(HMLE_df)
# divide the RPM values by the length of the gene in kilobase
# load in the conversion table:
convert<-read.table('~/Downloads/genesym_to_ensembl_genelengths.txt', header = F, sep = '\t')
head(convert)
# dimensions: 45471;3
dim(convert)
convert_HTE<- merge(data.frame('ensembl'=rownames(HMLE_df), HMLE_df), convert, by.x='ensembl', by.y='V1')
head(convert_HTE)
convert_HTE<-convert_HTE[!duplicated(convert_HTE$ensembl),]
# dimensions: 45471;8
head(convert_HTE)
dim(convert_HTE)
final_df<-convert_HTE[,c(2:7)]
head(final_df)
# normalize by gene length
for (l in 1:nrow(final_df)){
  final_df[l,]<-(as.numeric(convert_HTE[l,2:7])*1000)/convert_HTE[l,9]
}
rownames(final_df)<-convert_HTE$ensembl
head(final_df)

#  cat hallmark_interferon_alpha_response_gene_set.txt hallmark_interferon_gamma_response_gene_set.txt A.J.Scadden_nature_struct_bio.interferon_I_list|awk '!seen[$0]++' > ISG.txt
ISG<-read.table('~/Downloads/ISG.txt')
head(ISG)
combat_df<-data.frame('hgnc_symbol'=rownames(final_df), final_df)
#combat_df<-data.frame('hgnc_symbol'=symbol_HMLE$hgnc_symbol, symbol_HMLE[,-c(1:3)])
hnM_ISG<-merge(combat_df, ISG, by.x='hgnc_symbol', by.y='V1')
hnM_ISG<-hnM_ISG[!duplicated(hnM_ISG$hgnc_symbol),]
head(hnM_ISG)
# 195
dim(hnM_ISG)
ISG_input<-hnM_ISG[,-1]
rownames(ISG_input)<-hnM_ISG$hgnc_symbol
dim(ISG_input)
head(ISG_input)

################################################################################
for (i in 1:ncol(ISG_input)){
  ISG_input[,i]<-as.numeric(ISG_input[,i])
}
summary(ISG_input)
FC_input<-ISG_input
# 195; 9
dim(FC_input)
for (k in 1:nrow(FC_input)){
  FC_input[k,1]<-(ISG_input[k,1]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  FC_input[k,2]<-(ISG_input[k,2]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  FC_input[k,3]<-(ISG_input[k,3]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  FC_input[k,4]<-(ISG_input[k,4]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  FC_input[k,5]<-(ISG_input[k,5]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  FC_input[k,6]<-(ISG_input[k,6]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
}
mean_KD<-c()
for (i in 1:nrow(FC_input)){
  mean_KD<-c(mean_KD, mean(as.numeric(FC_input[i,])))
}
FC_input$mean_KD<-mean_KD
head(FC_input)
summary(FC_input)
new_rev_trim<-FC_input[rev(order(FC_input$mean_KD)),-ncol(FC_input)]
new_rev_trim<-na.omit(new_rev_trim)
new_rev_trim[new_rev_trim>=5]<-5
summary(new_rev_trim)
head(new_rev_trim)
rev_rev_trim<-new_rev_trim[rev(order(rowSums(new_rev_trim))),]
### heatmap visualization with viral-related interferon response genes highlighted
library(pheatmap)
library(RColorBrewer)
#install.packages('viridis')
library(viridis)
genelist<-c("MX1","IFI27","IFIT1","IFIT2","IFI3","DDX58","IFIH1","TLR3","CXCL10","IFI6",
            "OAS3","IRF7","ISG15","BST2")
labels<-rownames(rev_rev_trim)
labels[!labels %in% genelist] <- ""

makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}

cols <- makeColorRampPalette(c("royalblue", "white",    # distances 0 to 3 colored from white to red
                               "white", "red3"), # distances 3 to max(distmat) colored from green to black
                             0.17,
                             100)

pheatmap(mat=rev_rev_trim[1:100,],
         labels_row = labels[1:100],
         scale='none',
         color=cols,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         #color=colorRampPalette(c("royalblue","white","red","red3"))(n=30),
         cluster_rows = FALSE,
         cluster_cols = TRUE)
