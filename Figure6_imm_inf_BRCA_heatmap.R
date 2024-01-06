#############
#################################################################################
### create BRCA lcpmtmm matrix:
cancer_types = c('BRCA')
# define working directory
base_dir = '~/Downloads/'

# counts table
#counts_filename = sapply(cancer_types, function(x) paste(base_dir, 'GDC_', x,'_counts/',x,'_out_counts.txt',sep = ''))      
counts_filename = sapply(cancer_types, function(x) paste(base_dir,x,'_out_counts.txt',sep = ''))   
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
# #############################
# # # ###########################
# # # to extract normal samples:
# metadata<-read.csv('~/Downloads/clinical_data_BRCA.txt', sep = '\t')
# # PAM50Call_RNAseq
# # sample_type
# unique(metadata$PAM50Call_RNAseq)
# #LumA<-metadata[metadata$sample_type=="Solid Tissue Normal",c('sampleID','sample_type')]
# #LumA_ori<-metadata[metadata$PR_Status_nature2012=="Negative" & metadata$HER2_Final_Status_nature2012=="Negative" & metadata$ER_Status_nature2012=="Negative", c('sampleID','sample_type')]
# LumA_ori<-metadata[metadata$PAM50Call_RNAseq=="Basal",c('sampleID','sample_type')]
# LumA<-LumA_ori[LumA_ori$sample_type=="Primary Tumor",]
# # get rid of a metastatic sample
# dim(LumA)

#ACC_meta_new<-merge(ACC_meta, LumA, by.x="Sample.ID", by.y="sampleID")
ACC_meta_new<-ACC_meta
ACC_meta_new<-ACC_meta_new[!duplicated(ACC_meta_new$Sample.ID),]
dim(ACC_meta_new)
head(ACC_meta_new)
################################################################################
######lcpmtmm normalization
ACC_count<-read.table(counts_filename[[i]], sep = '\t', header = T)
ACC_count$Ensembl_gene_id<-substring(ACC_count$Ensembl_gene_id,1,15)
ACC_count<-ACC_count[!duplicated(ACC_count$Ensembl_gene_id),]
rownames(ACC_count)<-ACC_count$Ensembl_gene_id
colnames(ACC_count)<-gsub("X","", colnames(ACC_count))
#### get rid of low abundance genes
colnames(ACC_count)<-c('Ensembl_gene_id', gsub(".htseq.counts","",colnames(ACC_count[,-1])))
colnames(ACC_count)<-c('Ensembl_gene_id', gsub("\\.","-",colnames(ACC_count[,-1])))
dim(ACC_count)
head(ACC_count[1:3,1:3])
###########=====================================================================
# ACC_meta$File.Name<-gsub('\\.', '-', ACC_meta$File.Name)
# ACC_count_trim<-ACC_count[,c(ACC_meta$File.Name)]
ACC_meta_new$File.Name<-gsub('\\.', '-', ACC_meta_new$File.Name)
ACC_count_trim<-ACC_count[,c(ACC_meta_new$File.Name)]
head(ACC_count_trim[,1:5])
# 60483 139
dim(ACC_count_trim)
###
#BiocManager::install('DESeq2')
library('DESeq2')
#BiocManager::install('edgeR')
library('edgeR')

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
# 18926 139
dim(dds_kpt)
##################
ym<-DGEList(counts(dds_kpt))
ym<-calcNormFactors(ym, method="TMM")
lcpmtmm.rna<-cpm(ym, log=T, normalized.lib.sizes=T)
colnames(lcpmtmm.rna)<-gsub("X","",colnames(lcpmtmm.rna))
head(lcpmtmm.rna[,1:5])
###convert to gene symbol
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),
  values=rownames(lcpmtmm.rna),
  mart=mart)
head(genes)
dim(genes)
symbol<-merge(data.frame('ensembl_gene_id'=rownames(lcpmtmm.rna), lcpmtmm.rna), genes, by='ensembl_gene_id')
symbol<-symbol[!duplicated(symbol$hgnc_symbol),]
head(symbol[,c(1,(ncol(symbol)-1):ncol(symbol))])
dim(symbol)
output<-symbol[,-c(1,(ncol(symbol)-1):ncol(symbol))]
rownames(output)<-symbol$hgnc_symbol
write.table(output, '~/Downloads/BRCA_Basal_lcpmtmm.txt', sep = '\t', quote = F)

##############
library("htmltools")
#BiocManager::install('GSVA')
library("GSVA")
# output is empty
myoutf <- "~/Downloads/BRCA_Basal_tot_pathway_activity.xls"

filename = "~/Downloads/tot_imm_inf.gmt.txt"
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
#Set up the data set for running
# gene matrixs (colnames: cell ID)
# input file 
myinf <- "~/Downloads/BRCA_Basal_lcpmtmm.txt"
res <- read.table(myinf,sep="\t")
head(res[,1:4])
dim(res)
#Make the gene names identical
gene <- row.names(res)
#######################
# gene_CTC_EPI
for(i in 1 : length(all_pathway))
{
  cat("\r",i)
  com <- intersect(as.vector(unlist(all_pathway[[i]])), gene)
  all_pathway[[i]] <- com
}

gsva_es <- gsva(as.matrix(res), gset.idx.list=all_pathway, method="ssgsea",parallel.sz=4)
write.table(gsva_es,file = myoutf,row.names = TRUE,col.names = TRUE,sep = "\t",quote=F)
###
basal_df<-read.table('~/Downloads/BRCA_tot_hnrnpm_gsva.txt', sep = '\t', header = T)
basal_df
head(basal_df[,1:5])
dim(basal_df)
t_basal_df<-t(basal_df)
head(t_basal_df)
write.table(t_basal_df,'~/Downloads/tmp_basal_df.txt', sep = '\t', quote = F)
basal_order<-read.table('~/Downloads/tmp_basal_df.txt', sep = '\t', header = T)
basal_order$Cell_Types<-gsub('\\.','-',gsub('X','',basal_order$Cell_Types))
head(basal_order)

### 
convert_df<-read.table('~/Downloads/BRCA_Basal_heatmap/gdc_BRCA_counts.tsv', sep = '\t', header = T)
convert_df$File.Name<-substring(convert_df$File.Name, 1, 36)
head(convert_df)

basal_convert<-merge(basal_order, convert_df, by.x='Cell_Types', by.y='File.Name')
basal_convert<-basal_convert[!duplicated(basal_convert$Cell_Types),c(1:11,17)]
head(basal_convert)
basal_convert$Sample.ID<-substring(basal_convert$Sample.ID, 1, 12)

head(rs_df)
################################################################################
################################################################################
################################################################################
#library(TCGAbiolinks)
#install.packages('reshape2')
library(reshape2)
library(dplyr)
library(data.table)

library(survival)

#library(ranger)
#library(ggfortify)
library(survminer)
# ggfortify, ranger and survminer required but are not installed 


res_name = 'topcordsrnas_editing_survival' # replace


tcga_cdr_filename = '~/Downloads/TCGA-CDR_primary_sheet.txt'
# CDR can be downloaded from https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018
#tcga_cdr_filename = '~/Downloads/TCGA-CDR-SupplementalTableS1.xlsx'
# define working directory (work_dir)


###########################################################################################################

# function to trim both leading and trailing whitespace
trim = function(x) gsub("^\\s+|\\s+$", "", x)

make_survival_df = function(edit_df, clin_df){
  survival_edit_df = merge(edit_df, clin_df, by.x='sample', by.y='symbol')
  survival_edit_df$OS.time = as.numeric(survival_edit_df$OS.time)
  survival_edit_df$PFI.time = as.numeric(survival_edit_df$PFI.time)
  survival_edit_df$OS = as.numeric(survival_edit_df$OS)
  survival_edit_df$PFI = as.numeric(survival_edit_df$PFI)
  return(survival_edit_df)
}

test_survival_uni_cox_editing = function(survival_df, var_name){
  
  cox_os_fit = coxph(Surv(OS.time, OS) ~ get(var_name), 
                     data=survival_df[which(!is.na(survival_df[, c(get(var_name))])),])
  cox_os_ph_test = cox.zph(cox_os_fit)
  
  cox_pfi_fit = coxph(Surv(PFI.time, PFI) ~ get(var_name),
                      data=survival_df[which(!is.na(survival_df[, c(get(var_name))])),])
  cox_pfi_ph_test = cox.zph(cox_pfi_fit)
  
  #cat('Results from Cox PH OS:\n')
  #print(summary(cox_os_fit))
  cox_os_beta = summary(cox_os_fit)$coefficients[1, 1]
  cox_os_hr = summary(cox_os_fit)$coefficients[1, 2]
  cox_os_p = summary(cox_os_fit)$coefficients[1, 5]
  cat(paste('Cox_OS_pvalue ', cox_os_p, '\n', sep=''))
  #cat('Does the PH assumption hold for OS?\n')
  #print(cox_os_ph_test)
  
  #cat('Results from Cox PH PFI:\n')
  #print(summary(cox_pfi_fit))
  cox_pfi_beta = summary(cox_pfi_fit)$coefficients[1,1]
  cox_pfi_hr = summary(cox_pfi_fit)$coefficients[1,2]
  cox_pfi_p = summary(cox_pfi_fit)$coefficients[1, 5]
  cat(paste('Cox_PFI_pvalue ', cox_pfi_p, '\n', sep=''))
  #cat('Does the PH assumption hold for PFI?\n')
  #print(cox_pfi_ph_test)
  
  stat_df = data.frame(variable = rep(var_name, 2), endpoint = c('OS', 'PFI'), cox_pval = c(cox_os_p, cox_pfi_p),
                       cox_beta = c(cox_os_beta, cox_pfi_beta), cox_hr = c(cox_os_hr, cox_pfi_hr),
                       ph_pval = c(cox_os_ph_test$table[1,3], cox_pfi_ph_test$table[1,3]))
  return(list(cox_os_fit, cox_pfi_fit, stat_df, cox_os_ph_test, cox_pfi_ph_test))
}

test_survival_uni_cox_editing_string = function(survival_df, var_name){
  
  cox_os_fit = coxph(as.formula(paste('Surv(OS.time, OS) ~ ', var_name, sep='')), 
                     data=survival_df[which(!is.na(survival_df[, c(var_name)])),])
  cox_os_ph_test = cox.zph(cox_os_fit)
  
  cox_pfi_fit = coxph(as.formula(paste('Surv(PFI.time, PFI) ~ ', var_name, sep='')),
                      data=survival_df[which(!is.na(survival_df[, c(var_name)])),])
  cox_pfi_ph_test = cox.zph(cox_pfi_fit)
  
  #cat('Results from Cox PH OS:\n')
  #print(summary(cox_os_fit))
  cox_os_beta = summary(cox_os_fit)$coefficients[1, 1]
  cox_os_hr = summary(cox_os_fit)$coefficients[1, 2]
  cox_os_p = summary(cox_os_fit)$coefficients[1, 5]
  cat(paste('Cox_OS_pvalue ', cox_os_p, '\n', sep=''))
  #cat('Does the PH assumption hold for OS?\n')
  #print(cox_os_ph_test)
  
  #cat('Results from Cox PH PFI:\n')
  #print(summary(cox_pfi_fit))
  cox_pfi_beta = summary(cox_pfi_fit)$coefficients[1,1]
  cox_pfi_hr = summary(cox_pfi_fit)$coefficients[1,2]
  cox_pfi_p = summary(cox_pfi_fit)$coefficients[1, 5]
  cat(paste('Cox_PFI_pvalue ', cox_pfi_p, '\n', sep=''))
  #cat('Does the PH assumption hold for PFI?\n')
  #print(cox_pfi_ph_test)
  
  stat_df = data.frame(variable = rep(var_name, 2), endpoint = c('OS', 'PFI'), cox_pval = c(cox_os_p, cox_pfi_p),
                       cox_beta = c(cox_os_beta, cox_pfi_beta), cox_hr = c(cox_os_hr, cox_pfi_hr),
                       ph_pval = c(cox_os_ph_test$table[1,3], cox_pfi_ph_test$table[1,3]))
  return(list(cox_os_fit, cox_pfi_fit, stat_df, cox_os_ph_test, cox_pfi_ph_test))
}


###########################################################################################################

# prepare inputs

library(dplyr)
tcga_cdr_df = read.table(tcga_cdr_filename, sep='\t', header=T, stringsAsFactors = F, check.names = F, comment.char='') #, fill=T)
tcga_cdr_df = tcga_cdr_df %>% dplyr::filter(Redaction != "Redacted")

# paste <(awk '{print $1}' BRCA_MED15_Apr/CryEx_MED15.txt) <(awk '{print $4}' BRCA_MED15_Apr/CryEx_MED15.txt) <(awk '{print $4}' BRCA_RBM34_Apr/CryEx_RBM34.txt) <(awk '{print $5}' BRCA_TRAPPC10_Apr/CryEx_TRAPPC10.txt) <(awk '{print $5}' BRCA_LRP11_Apr/CryEx_LRP11.txt) > CryEx_tot4.txt
lq_TRAPPC10_Basal_tmp<-read.table('~/Downloads/BRCA_Basal_heatmap/CryEx_tot4.txt', sep = '\t', header = F)
lq_TRAPPC10_Basal_tmp$V1<-substring(lq_TRAPPC10_Basal_tmp$V1, 6, 41)
head(lq_TRAPPC10_Basal_tmp)
colnames(lq_TRAPPC10_Basal_tmp)<-c('File.ID','MED15','RBM34','TRAPPC10','LRP11')

converter<-read.table('~/Downloads/BRCA_Basal_heatmap/gdc_BRCA_BAM_Apr.tsv', sep = '\t', header = T)
head(converter)


input_BRCA<-merge(lq_TRAPPC10_Basal_tmp, converter, by='File.ID')
input_BRCA<-input_BRCA[!duplicated(input_BRCA$File.ID),]
head(input_BRCA)
dim(input_BRCA)
input_BRCA_res<-merge(input_BRCA, res_convert, by='Case.ID')
input_BRCA_res<-input_BRCA_res[!duplicated(input_BRCA_res$Case.ID),]

#write.table(input_BRCA, '~/Downloads/input_BRCA.txt', sep = '\t', row.names = F, quote = F)




# lq_TRAPPC10_Basal<-read.table('~/Downloads/ISGs_PSIs_hnM_BRCA_Basal.txt', sep = '\t', header = T)
# lq_TRAPPC10_Basal$Sample.ID<-substring(lq_TRAPPC10_Basal$Sample.ID, 1, 12)
# head(lq_TRAPPC10_Basal)
# dim(lq_TRAPPC10_Basal)
# l_hnM<-lq_TRAPPC10_Basal[lq_TRAPPC10_Basal$hnM< mean(lq_TRAPPC10_Basal$hnM),]
# h_hnM<-lq_TRAPPC10_Basal[lq_TRAPPC10_Basal$hnM> mean(lq_TRAPPC10_Basal$hnM),]
# label_hnM<-c(rep('low',nrow(l_hnM)), rep('high',nrow(h_hnM)))
# input<-data.frame(rbind(l_hnM, h_hnM))
# input$label_hnM=label_hnM


survival_edit_df = merge(input_BRCA_res, tcga_cdr_df, by.x='Case.ID', by.y='bcr_patient_barcode')
dim(survival_edit_df)
summary(survival_edit_df)
survival_edit_df<-na.omit(survival_edit_df)
survival_edit_df$OS.time = as.numeric(survival_edit_df$OS.time)
survival_edit_df$PFI.time = as.numeric(survival_edit_df$PFI.time)
survival_edit_df$OS = as.numeric(survival_edit_df$OS)
survival_edit_df$PFI = as.numeric(survival_edit_df$PFI)




# ##############################################################################################################################################
##############################################################################################################################################
# Multivariate Cox regression


filt_survival_df = survival_edit_df # replace editing_mean with desired feature
#filt_survival_df$race = sapply(filt_survival_df$race, function(x) ifelse(x == "[Not Available]", NA,
#    ifelse(x == "[Not Evaluated]", NA, 
#     ifelse(x == "[Unknown]", NA, x))))
#filt_survival_df$age = as.numeric(filt_survival_df$age_at_initial_pathologic_diagnosis)
filt_survival_df<-na.omit(filt_survival_df)
head(filt_survival_df)
dim(filt_survival_df)

#TRAPPC10 + LRP11 + 
library('survival')

risk_score_lou<-c()
for (i in 1:nrow(filt_survival_df)){
  filt_survival_df_trim<-filt_survival_df[-i,]
  mult_cox_os_fit = coxph(Surv(OS.time, OS) ~ TRAPPC10 + LRP11 + RBM34 + MED15 + hnM,   # replace editing_mean with desired feature (and covariates as desired)
                          data=filt_survival_df_trim)
  risk_score<-predict(mult_cox_os_fit, filt_survival_df, type = "risk")
  risk_score_lou<-c(risk_score_lou, risk_score[i])
}


# forest plot
filt_survival_df<-na.omit(filt_survival_df)
dim(filt_survival_df)
#p = ggforest(mult_cox_os_fit, fontsize = 1, main = 'Hazard Ratios for OS', 
#             data=filt_survival_df) #, cpositions = c(0.02, 0.15, 0.2)) ggtheme = theme_survminer(),
#print(p)



filt_survival_df$risk_score<-risk_score_lou
rs_low<-filt_survival_df[filt_survival_df$risk_score< mean(filt_survival_df$risk_score),]
rs_high<-filt_survival_df[filt_survival_df$risk_score> mean(filt_survival_df$risk_score),]
# as.numeric(quantile(as.numeric(lcpmtmm.rna[rownames(lcpmtmm.rna)=="ENSG00000099783",]), probs = c(0.25, 0.75))[1])

label_rs<-c(rep('low',nrow(rs_low)), rep('high',nrow(rs_high)))
rs_df<-data.frame(rbind(rs_low, rs_high))
rs_df$label_rs<-label_rs
head(rs_df)

fit_OS<-survfit(Surv(OS.time, OS) ~ label_rs, data = rs_df)

fit_PFI<-survfit(Surv(PFI.time, PFI) ~ label_rs, data = rs_df)
#install.packages('survminer')
#library('survminer')
pp_OS=ggsurvplot(fit_OS, data = rs_df, pval=TRUE)
pp_PFI=ggsurvplot(fit_PFI, data = rs_df, pval=TRUE)
pp_OS
pp_PFI

###
basal_df<-read.table('~/Downloads/BRCA_tot_hnrnpm_gsva.txt', sep = '\t', header = T)
basal_df
head(basal_df[,1:5])
dim(basal_df)
t_basal_df<-t(basal_df)
head(t_basal_df)
write.table(t_basal_df,'~/Downloads/tmp_basal_df.txt', sep = '\t', quote = F)
basal_order<-read.table('~/Downloads/tmp_basal_df.txt', sep = '\t', header = T)
basal_order$Cell_Types<-gsub('\\.','-',gsub('X','',basal_order$Cell_Types))
head(basal_order)

converter_ori<-read.table('~/Downloads/BRCA_Basal_heatmap/gdc_BRCA_counts.tsv', sep = '\t', header = T)
converter_ori$File.Name<-substring(converter_ori$File.Name, 1, 36)
head(converter_ori)
dim(converter_ori)

basal_converter<-merge(basal_order, converter_ori, by.x='Cell_Types', by.y='File.Name')
basal_converter<-basal_converter[!duplicated(basal_converter$Cell_Types),]
head(basal_converter)
dim(basal_converter)
write.table(basal_converter, '~/Downloads/basal_conv.txt', sep = '\t', quote = F, row.names = F)

### 

head(rs_df)

basal_rs<-merge(basal_converter, rs_df, by='Case.ID')
# basal_rs<-basal_rs[!duplicated(basal_rs$Case.ID),c(2:12,18:22,56)]
# dim(basal_rs)
head(basal_rs)

basal_rs_order<-basal_rs[rev(order(basal_rs$MED15)),]
head(basal_rs_order)
dim(basal_rs_order)

t_basal_rs<-t(basal_rs_order)
head(t_basal_rs[,1:5])

write.table(t_basal_rs, '~/Downloads/t_basal_rs.txt', sep = '\t', quote = F)
mm_input<-read.csv('~/Downloads/t_basal_rs.txt', sep = '\t', header = T)
#mm_input<-read.table('~/Downloads/t_df_order.txt', sep = '\t', header = T)
head(mm_input[,1:5])
rownames(mm_input)<-mm_input$Case.ID

mm_trim<-mm_input[-c(11:17, 19:nrow(mm_input)),-1]
for (i in 1:ncol(mm_trim)){
  mm_trim[,i]<-as.numeric(mm_trim[,i])
}
summary(mm_trim[,1:3])

library(pheatmap)
library(RColorBrewer)
#install.packages('viridis')
library(viridis)


makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}

cols <- makeColorRampPalette(c("blue", "white",    # distances 0 to 3 colored from white to red
                               "white", "red"), # distances 3 to max(distmat) colored from green to black
                             0.5,
                             100)

pheatmap(mm_trim, 
         scale = 'row',
         color=cols,
         cluster_rows = F,
         cluster_cols = F
)

mm_trim_z<-mm_trim
for (l in 1:nrow(mm_trim)){
  for (k in 1:ncol(mm_trim)){
    mm_trim_z[l,k]<-(mm_trim[l,k]-mean(as.numeric(mm_trim[l,])))/sd(as.numeric(mm_trim[l,]))
  }
}
pheatmap(mm_trim_z,
         color=cols,
         cluster_rows = F,
         cluster_cols = F)

mm_trim_z_z<-mm_trim_z
mm_trim_z_z[mm_trim_z_z>= 2]<-2
mm_trim_z_z[mm_trim_z_z<= -1]<- -1
mm_input_z<-data.frame(rbind(mm_trim_z_z, mm_trim[11,]))
mm_input_z_z<-mm_input_z
tail(mm_input_z_z[,1:5])
for (z in 1:ncol(mm_input_z_z)){
  mm_input_z_z[11,z]<-(mm_input_z[11,z]-mean(as.numeric(mm_input_z[11,])))/sd(as.numeric(mm_input_z[11,]))
}
tail(mm_input_z_z[,1:5])
summary(as.numeric(mm_input_z_z[11,]))

mm_input_z_z_z<-mm_input_z_z
mm_input_z_z_z[mm_input_z_z>= 2]<-2
mm_input_z_z_z[mm_input_z_z<= -1]<- -1
pheatmap(mm_input_z_z_z[-c(4:5,7),],
         color=cols,
         cluster_rows = F,
         cluster_cols = F)
pheatmap(mm_input_z_z_z[-c(4:5,7),1:274],
         color=cols,
         cluster_rows = F,
         cluster_cols = F)
pheatmap(mm_input_z_z_z[-c(4:5,7),820:1094],
         color=cols,
         cluster_rows = F,
         cluster_cols = F)

################################################################################
t_mm_input<-t(mm_input)
write.table(t_mm_input, '~/Downloads/t_mm_input.txt', sep = '\t', quote = F, row.names = F)
check_mm<-read.csv('~/Downloads/t_mm_input.txt', sep = '\t')
check_mm$Sample.ID.x<-substring(check_mm$Sample.ID.x, 1, 12)
head(check_mm[,1:3])
dim(check_mm)

#################################################################################################
res_trim<-data.frame('samples'=colnames(res), 'hnM'=as.numeric(res[rownames(res)=="HNRNPM",]))
res_trim$samples<-gsub('\\.','-',gsub('X','', res_trim$samples))
head(res_trim)
head(convert_df)

res_convert<-merge(res_trim, convert_df, by.x='samples', by.y='File.Name')
res_convert<-res_convert[!duplicated(res_convert$samples),]
head(res_convert)


### 

wilcox.test(as.numeric(check_mm[check_mm$label_rs=="low",2]), as.numeric(check_mm[check_mm$label_rs=="high",2]))


wilcox.test(as.numeric(check_mm[check_mm$MED15<=mean(as.numeric(check_mm$MED15)),2]), as.numeric(check_mm[check_mm$MED15>=mean(as.numeric(check_mm$MED15)),2]))

t.test(as.numeric(check_mm[check_mm$MED15<=mean(as.numeric(check_mm$MED15)),2]), as.numeric(check_mm[check_mm$MED15>=mean(as.numeric(check_mm$MED15)),2]))


wilcox.test(as.numeric(check_mm[check_mm$RBM34<=mean(as.numeric(check_mm$RBM34)),2]), as.numeric(check_mm[check_mm$RBM34>=mean(as.numeric(check_mm$RBM34)),2]))

t.test(as.numeric(check_mm[check_mm$RBM34<=mean(as.numeric(check_mm$RBM34)),2]), as.numeric(check_mm[check_mm$RBM34>=mean(as.numeric(check_mm$RBM34)),2]))


mean(as.numeric(check_mm[check_mm$label_rs=="low",2]))
mean(as.numeric(check_mm[check_mm$label_rs=="high",2]))


mean(as.numeric(check_mm[check_mm$RBM34<=mean(as.numeric(check_mm$RBM34)),4]))
mean(as.numeric(check_mm[check_mm$RBM34>=mean(as.numeric(check_mm$RBM34)),4]))
##########################
### DLBC
#################################################################################
cancer_types = c('TGCT')
# define working directory
base_dir = '~/Downloads/'

# counts table
#counts_filename = sapply(cancer_types, function(x) paste(base_dir, 'GDC_', x,'_counts/',x,'_out_counts.txt',sep = ''))      
counts_filename = sapply(cancer_types, function(x) paste(base_dir,x,'_out_counts.txt',sep = ''))   
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
################################################################################
######lcpmtmm normalization
ACC_count<-read.table(counts_filename[[i]], sep = '\t', header = T)
ACC_count$Ensembl_gene_id<-substring(ACC_count$Ensembl_gene_id,1,15)
ACC_count<-ACC_count[!duplicated(ACC_count$Ensembl_gene_id),]
rownames(ACC_count)<-ACC_count$Ensembl_gene_id
colnames(ACC_count)<-gsub("X","", colnames(ACC_count))
#### get rid of low abundance genes
colnames(ACC_count)<-c('Ensembl_gene_id', gsub(".htseq.counts","",colnames(ACC_count[,-1])))
colnames(ACC_count)<-c('Ensembl_gene_id', gsub("\\.","-",colnames(ACC_count[,-1])))
dim(ACC_count)
head(ACC_count[1:3,1:3])
###########=====================================================================
ACC_meta$File.Name<-gsub('\\.', '-', ACC_meta$File.Name)
ACC_count_trim<-ACC_count[,c(ACC_meta$File.Name)]
head(ACC_count_trim[,1:5])
# 60483 139
dim(ACC_count_trim)
###
#BiocManager::install('DESeq2')
#library('DESeq2')
#BiocManager::install('edgeR')
#library('edgeR')

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
# 18926 139
dim(dds_kpt)
##################
ym<-DGEList(counts(dds_kpt))
ym<-calcNormFactors(ym, method="TMM")
lcpmtmm.rna<-cpm(ym, log=T, normalized.lib.sizes=T)
colnames(lcpmtmm.rna)<-gsub("X","",colnames(lcpmtmm.rna))
head(lcpmtmm.rna[,1:5])

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
# 46
length(lq)
# 46
length(hq)
################################################################################
count_lq_hq<-counts(dds_kpt)[,c(lq,hq)]
head(count_lq_hq)
lq_hq<-data.frame(cbind('sample'=colnames(count_lq_hq), 'label'=c(rep('lq',length(lq)), rep('hq',length(hq)))))
head(lq_hq)

converter<-read.table('~/Downloads/UVM_dsRNA/gdc_UVM_counts.tsv', sep = '\t', header = T)
converter$File.Name<-substring(converter$File.Name, 1, 36)
head(converter)

lq_hq_converter<-merge(lq_hq, converter, by.x='sample', by.y='File.Name')
head(lq_hq_converter)

intermediate_label<-lq_hq_converter[,c('label','Case.ID')]
head(intermediate_label)

CryEx_LRP11<-read.table('~/Downloads/UVM_dsRNA/CryEx_MED15.txt', sep = '\t')
CryEx_LRP11$V1<-substring(CryEx_LRP11$V1,5,40)
CryEx_LRP11<-na.omit(CryEx_LRP11)
head(CryEx_LRP11)

converter_Apr<-read.table('~/Downloads/UVM_dsRNA/gdc_UVM_BAM_Apr.tsv', sep = '\t', header = T)
head(converter_Apr)

CryEx_LRP11_converter<-merge(CryEx_LRP11, converter_Apr, by.x='V1', by.y='File.ID')
CryEx_LRP11_converter<-CryEx_LRP11_converter[!duplicated(CryEx_LRP11_converter$V1),]
head(CryEx_LRP11_converter)

intermediate_LRP11<-CryEx_LRP11_converter[,c('Case.ID','V4')]
head(intermediate_LRP11)

lq_hq_LRP11<-merge(intermediate_label, intermediate_LRP11, by='Case.ID')
lq_hq_LRP11<-lq_hq_LRP11[!duplicated(lq_hq_LRP11$Case.ID),]
head(lq_hq_LRP11)

#boxplot(lq_hq_LRP11[lq_hq_LRP11$label=="lq",]$V5, lq_hq_LRP11[lq_hq_LRP11$label=="hq",]$V5)

boxplot(log2(lq_hq_LRP11[lq_hq_LRP11$label=="lq",]$V4+0.001), log2(lq_hq_LRP11[lq_hq_LRP11$label=="hq",]$V4+0.001))








