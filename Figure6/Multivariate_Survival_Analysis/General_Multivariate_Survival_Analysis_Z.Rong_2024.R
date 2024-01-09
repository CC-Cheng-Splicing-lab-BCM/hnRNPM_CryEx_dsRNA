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

library(dplyr)
tcga_cdr_df = read.table(tcga_cdr_filename, sep='\t', header=T, stringsAsFactors = F, check.names = F, comment.char='') #, fill=T)
tcga_cdr_df = tcga_cdr_df %>% dplyr::filter(Redaction != "Redacted")


# ##############################################################################################################################################
##############################################################################################################################################
# Multivariate Cox regression


filt_survival_df = survival_edit_df # replace editing_mean with desired feature
filt_survival_df<-na.omit(filt_survival_df)
head(filt_survival_df)
dim(filt_survival_df)

library('survival')

risk_score_lou<-c()
for (i in 1:nrow(filt_survival_df)){
  filt_survival_df_trim<-filt_survival_df[-i,]
  # TRAPPC10, LRP11, RBM34 and MED15 are the CryEx PSIs calculated from TCGA cancer samples 
  # hnM is the normalized hnRNPM expression level calculated from TCGA cancer samples
  mult_cox_os_fit = coxph(Surv(OS.time, OS) ~ TRAPPC10 + LRP11 + RBM34 + MED15 + hnM,   # replace editing_mean with desired feature (and covariates as desired)
                          data=filt_survival_df_trim)
  # Risk scores is calculated
  risk_score<-predict(mult_cox_os_fit, filt_survival_df, type = "risk")
  risk_score_lou<-c(risk_score_lou, risk_score[i])
}


# forest plot
filt_survival_df<-na.omit(filt_survival_df)
dim(filt_survival_df)

filt_survival_df$risk_score<-risk_score_lou
rs_low<-filt_survival_df[filt_survival_df$risk_score< mean(filt_survival_df$risk_score),]
rs_high<-filt_survival_df[filt_survival_df$risk_score> mean(filt_survival_df$risk_score),]

label_rs<-c(rep('low',nrow(rs_low)), rep('high',nrow(rs_high)))
rs_df<-data.frame(rbind(rs_low, rs_high))
rs_df$label_rs<-label_rs
head(rs_df)

fit_OS<-survfit(Surv(OS.time, OS) ~ label_rs, data = rs_df)

fit_PFI<-survfit(Surv(PFI.time, PFI) ~ label_rs, data = rs_df)
#install.packages('survminer')
library('survminer')
pp_OS=ggsurvplot(fit_OS, data = rs_df, pval=TRUE)
pp_PFI=ggsurvplot(fit_PFI, data = rs_df, pval=TRUE)
pp_OS
pp_PFI