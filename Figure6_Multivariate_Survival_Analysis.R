# input file 
myinf <- "~/Downloads/BRCA_Basal_lcpmtmm.txt"
res <- read.table(myinf,sep="\t")
head(res[,1:4])
dim(res)
#Make the gene names identical
gene <- row.names(res)

### 
convert_df<-read.table('~/Downloads/BRCA_Basal_heatmap/gdc_BRCA_counts.tsv', sep = '\t', header = T)
convert_df$File.Name<-substring(convert_df$File.Name, 1, 36)
head(convert_df)

#################################################################################################
res_trim<-data.frame('samples'=colnames(res), 'hnM'=as.numeric(res[rownames(res)=="HNRNPM",]))
res_trim$samples<-gsub('\\.','-',gsub('X','', res_trim$samples))
head(res_trim)
head(convert_df)

res_convert<-merge(res_trim, convert_df, by.x='samples', by.y='File.Name')
res_convert<-res_convert[!duplicated(res_convert$samples),]
head(res_convert)

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



library(dplyr)
tcga_cdr_df = read.table(tcga_cdr_filename, sep='\t', header=T, stringsAsFactors = F, check.names = F, comment.char='') #, fill=T)
tcga_cdr_df = tcga_cdr_df %>% dplyr::filter(Redaction != "Redacted")


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

survival_edit_df = merge(input_BRCA_res, tcga_cdr_df, by.x='Case.ID', by.y='bcr_patient_barcode')
dim(survival_edit_df)
summary(survival_edit_df)
survival_edit_df<-na.omit(survival_edit_df)
survival_edit_df$OS.time = as.numeric(survival_edit_df$OS.time)
survival_edit_df$PFI.time = as.numeric(survival_edit_df$PFI.time)
survival_edit_df$OS = as.numeric(survival_edit_df$OS)
survival_edit_df$PFI = as.numeric(survival_edit_df$PFI)

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
