################################################################################
### ADAR ### A-I editing sites ###
DE_AI<-read.table('~/Downloads/Total_summary.DE_sites_annotated.txt', sep = '\t', header = T)
head(DE_AI)
dim(DE_AI)
DE_AI_ECL_EKD<-DE_AI[DE_AI$Comparison=="ECL_vs_EKD" & DE_AI$DE.NS=="DE" & DE_AI$region=="intron",]
dim(DE_AI_ECL_EKD)
head(DE_AI_ECL_EKD)
write.table(DE_AI_ECL_EKD, '~/Downloads/DE_AI_ECL_EKD.txt', sep = '\t', quote = F)
### EKD
EKD<-read.table('~/Downloads/all_editing_sites.EKD.minCov_5.no_batch_mean.both_rep_cov.tab', sep = '\t', header = T)
head(EKD[,1:3])
dim(EKD)
ECL<-read.table('~/Downloads/all_editing_sites.ECL.minCov_5.no_batch_mean.both_rep_cov.tab', sep = '\t', header = T)
head(ECL[,1:3])
dim(ECL)
ECL_EKD<-merge(ECL, EKD, by='Site_ID', all = T)
ECL_EKD<-ECL_EKD[!duplicated(ECL_EKD$Site_ID),]
head(ECL_EKD)
dim(ECL_EKD)
ECL_EKD[is.na(ECL_EKD)]=0
summary(ECL_EKD)
library(tidyverse)
BiocManager::install('gapminder')
library(gapminder)
df = gapminder %>%
  filter(year %in% c(1952,2007)) %>%
  filter(continent %in% c("Americas")) %>%
  select(country,year,lifeExp)%>%
  mutate(paired = rep(1:(n()/2),each=2),
         year=factor(year))

df %>%
  ggplot(aes(year,lifeExp, fill=year)) +
  geom_boxplot() +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=year,group=paired),size=2,shape=21, position = position_dodge(0.2)) +
  theme(legend.position = "none")

plot(ECL_EKD$MeanRatio.x, ECL_EKD$MeanRatio.y)
write.table(ECL_EKD, '~/Downloads/ECL_EKD.txt', sep = '\t', quote = F, row.names = F)
# cat ECL_EKD.txt|awk '{print $1}'|cut -f1 -d ':' > chr
# cat ECL_EKD.txt|awk '{print $1}'|cut -f2 -d ':' > start
# paste chr start ECL_EKD.txt> add_ECL_EKD.txt
add_ECL_EKD<-read.table('~/Downloads/add_ECL_EKD.txt', sep = '\t', header = T)
head(add_ECL_EKD)
dim(add_ECL_EKD)
End<-c()
for (i in 1:nrow(add_ECL_EKD)){
  End<-c(End, (add_ECL_EKD$Start[i]+1))
}
add_ECL_EKD$End<-End
head(add_ECL_EKD)
write.table(add_ECL_EKD, '~/Downloads/add_ECL_EKD_add.bed', sep = '\t', quote = F, row.names = F)


LRP11_AluJb<-add_ECL_EKD[add_ECL_EKD$chr=="chr6" & add_ECL_EKD$Start>150178386 & add_ECL_EKD$Start<150178671,]$MeanRatio.x
LRP11_AluSz<-add_ECL_EKD[add_ECL_EKD$chr=="chr6" & add_ECL_EKD$Start>150180534 & add_ECL_EKD$Start<150180831,]$MeanRatio.x

LRP11_AluJb_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr6" & add_ECL_EKD$Start>150178386 & add_ECL_EKD$Start<150178671,]$MeanRatio.y
LRP11_AluSz_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr6" & add_ECL_EKD$Start>150180534 & add_ECL_EKD$Start<150180831,]$MeanRatio.y

boxplot(c(LRP11_AluJb, LRP11_AluSz), c(LRP11_AluSz_KD, LRP11_AluJb_KD))

colnames(add_ECL_EKD)<-c('chr','Start','Site_ID','MeanRatio.x','MeanCov.x','MeanRatio.y','MeanCov.y')
RBM34_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235305599 & add_ECL_EKD$Start<235305892,]
RBM34_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235305961 & add_ECL_EKD$Start<235306262,]
RBM34_3<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235306769 & add_ECL_EKD$Start<235307050,]
RBM34_4<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235307215 & add_ECL_EKD$Start<235307487,]
RBM34_5<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235307919 & add_ECL_EKD$Start<235308072,]
RBM34_6<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235308076 & add_ECL_EKD$Start<235308385,]
RBM34_7<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235308597 & add_ECL_EKD$Start<235308892,]
RBM34_8<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235309105 & add_ECL_EKD$Start<235309267,]
RBM34_9<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235309933 & add_ECL_EKD$Start<235310235,]
RBM34_10<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235309282 & add_ECL_EKD$Start<235309560,]
RBM34_11<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235310272 & add_ECL_EKD$Start<235310580,]
RBM34_ADAR<-data.frame(rbind(RBM34_1, RBM34_2, RBM34_3, RBM34_4, RBM34_5, RBM34_6, RBM34_7, 
                             RBM34_8, RBM34_9, RBM34_10, RBM34_11
                             ))
dim(RBM34_ADAR)
head(RBM34_ADAR)
# cat add_ECL_EKD_add.bed|awk '($1=="chr22")'|awk '($2>20921104)'|awk '($2<20936898)'|awk '{print $1"\t"$2"\t"($2+1)"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'> MED15_IGV_ADAR_real.bed
MED15_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr22" & add_ECL_EKD$Start>20927249 & add_ECL_EKD$Start<20927554,]
MED15_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr22" & add_ECL_EKD$Start>20926927 & add_ECL_EKD$Start<20927229,]
MED15_ADAR<-data.frame(rbind(MED15_1, MED15_2))
dim(MED15_ADAR)
head(MED15_ADAR)

# cat add_ECL_EKD_add.bed|awk '($1=="chr19")'|awk '($2>33183580)'|awk '($2<33200091)'|awk '{print $1"\t"$2"\t"($2+1)"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'> NUDT19_IGV_ADAR_real.bed
NUDT19_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33191503 & add_ECL_EKD$Start<33196790,]
NUDT19_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33191503 & add_ECL_EKD$Start<33196790,]
NUDT19_ADAR<-data.frame(rbind(NUDT19_1, NUDT19_2))
dim(NUDT19_ADAR)
head(NUDT19_ADAR)


TRAPPC10_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr21" & add_ECL_EKD$Start>45461255 & add_ECL_EKD$Start<45461551,]
TRAPPC10_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr21" & add_ECL_EKD$Start>45459466 & add_ECL_EKD$Start<45459761,]
TRAPPC10_ADAR<-data.frame(rbind(TRAPPC10_1, TRAPPC10_2
))
dim(TRAPPC10_ADAR)
head(TRAPPC10_ADAR)
#TRAPPC10_ADAR<-RBM34_ADAR
TRAPPC10_ADAR<-RBM34_ADAR
#TRAPPC10_ADAR<-NUDT19_ADAR
df_TRAPPC10<-data.frame('year'=c(rep('Ctrl',nrow(TRAPPC10_ADAR)), rep('KD',nrow(TRAPPC10_ADAR))),
                     'MeanRatio'=c(TRAPPC10_ADAR$MeanRatio.x, TRAPPC10_ADAR$MeanRatio.y),
                     'paired'=c(c(1:nrow(TRAPPC10_ADAR)), c(1:nrow(TRAPPC10_ADAR))))
head(df_RBM34)
dim(df_RBM34)
df_TRAPPC10 %>%
  ggplot(aes(year,MeanRatio, fill=year)) +
  geom_boxplot() +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=year,group=paired),size=2,shape=21, position = position_dodge(0.2)) +
  theme(legend.position = "none")
############################################
### Total:
df_tot<-data.frame('year'=c(rep('Ctrl', nrow(add_ECL_EKD)), rep('KD', nrow(add_ECL_EKD))),
                   'MeanRatio'=c(add_ECL_EKD$MeanRatio.x, add_ECL_EKD$MeanRatio.y),
                   'paired'=c(c(1:nrow(add_ECL_EKD)), c(1:nrow(add_ECL_EKD))))
df_tot %>%
  ggplot(aes(year,MeanRatio, fill=year)) +
  geom_boxplot() +
  #geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=year,group=paired),size=0.3,shape=21, position = position_dodge(0.2)) +
  theme(legend.position = "none")
###
# bedtools intersect -a CryEx_BED6_gencode.bed -b ECL_EKD.bed -wb|awk '!seen[$0]++'|awk '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'|awk '!seen[$0]++'> CryEx_ECL_EKD.bed
# bedtools intersect -a AS_BED6_gencode.bed -b ECL_EKD.bed -wb|awk '!seen[$0]++'|awk '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'|awk '!seen[$0]++'> AS_ECL_EKD.bed 
# bedtools intersect -a IR_BED6_gencode.bed -b ECL_EKD.bed -wb|awk '!seen[$0]++'|awk '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'|awk '!seen[$0]++'> IR_ECL_EKD.bed
# bedtools intersect -a bkgd_BED6_gencode.bed -b ECL_EKD.bed -wb|awk '!seen[$0]++'|awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}'|awk '!seen[$0]++'> bkgd_ECL_EKD.bed
# cat bkgd_ECL_EKD.bed IR_ECL_EKD.bed AS_ECL_EKD.bed CryEx_ECL_EKD.bed|awk '!seen[$0]++'> Exon_ECL_EKD.bed
# bedtools intersect -a ECL_EKD.bed -b Exon_ECL_EKD.bed -v|awk '!seen[$0]++'> nonExon_ECL_EKD.bed
library(dplyr)
library(ggplot2)
### Exon:
Exon_ECL_EKD<-read.table('~/Downloads/Exon_ECL_EKD.bed', sep = '\t')
head(Exon_ECL_EKD)
dim(Exon_ECL_EKD)

df_exon<-data.frame('year'=c(rep('Ctrl', nrow(Exon_ECL_EKD)), rep('KD', nrow(Exon_ECL_EKD))),
                   'MeanRatio'=c(Exon_ECL_EKD$V5, Exon_ECL_EKD$V7),
                   'paired'=c(c(1:nrow(Exon_ECL_EKD)), c(1:nrow(Exon_ECL_EKD))))
df_exon %>%
  ggplot(aes(year,MeanRatio, fill=year)) +
  geom_boxplot() +
  #geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=year,group=paired),size=0.3,shape=21, position = position_dodge(0.2)) +
  theme(legend.position = "none")
wilcox.test(Exon_ECL_EKD$V5, Exon_ECL_EKD$V7)
### non-Exon:
nonExon_ECL_EKD<-read.table('~/Downloads/nonExon_ECL_EKD.bed', sep = '\t')
head(nonExon_ECL_EKD)
dim(nonExon_ECL_EKD)

df_nonexon<-data.frame('year'=c(rep('Ctrl', nrow(nonExon_ECL_EKD)), rep('KD', nrow(nonExon_ECL_EKD))),
                    'MeanRatio'=c(nonExon_ECL_EKD$V5, nonExon_ECL_EKD$V7),
                    'paired'=c(c(1:nrow(nonExon_ECL_EKD)), c(1:nrow(nonExon_ECL_EKD))))
df_nonexon %>%
  ggplot(aes(year,MeanRatio, fill=year)) +
  geom_boxplot() +
  #geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=year,group=paired),size=0.3,shape=21, position = position_dodge(0.2)) +
  theme(legend.position = "none")
wilcox.test(nonExon_ECL_EKD$V5, nonExon_ECL_EKD$V7)
##################################################################
df_df<-data.frame('year'=c(rep('Ctrl_Exon',nrow(Exon_ECL_EKD)), rep('KD_Exon',nrow(Exon_ECL_EKD)), rep('Ctrl_nonExon',nrow(nonExon_ECL_EKD)),
                           rep('KD_nonExon',nrow(nonExon_ECL_EKD))),
                  'MeanRatio'=c(Exon_ECL_EKD$V5, Exon_ECL_EKD$V7, nonExon_ECL_EKD$V5, nonExon_ECL_EKD$V7),
                  'paired'=c(c(1:nrow(Exon_ECL_EKD)), c(1:nrow(Exon_ECL_EKD)), c((nrow(Exon_ECL_EKD)+1):(nrow(Exon_ECL_EKD)+nrow(nonExon_ECL_EKD))), 
                             c((nrow(Exon_ECL_EKD)+1):(nrow(Exon_ECL_EKD)+nrow(nonExon_ECL_EKD)))))
df_df %>%
  ggplot(aes(year,MeanRatio, fill=year))+
  geom_boxplot()+
  geom_point(aes(fill=year,group=paired), size=0.3,shape=21,position = position_dodge(0.2))+
  theme(legend.position = "none")
wilcox.test(Exon_ECL_EKD$V5,nonExon_ECL_EKD$V5)
###################################################################
### ADAR editing levels:
input_CryEx<-read.table('~/Downloads/CryEx_ECL_EKD.bed', sep = '\t', header = F)
head(input_CryEx)
dim(input_CryEx)
input_IR<-read.table('~/Downloads/IR_ECL_EKD.bed', sep = '\t', header = F)
head(input_IR)
dim(input_IR)
input_bkgd<-read.table('~/Downloads/bkgd_ECL_EKD.bed', sep = '\t', header = F)
head(input_bkgd)
dim(input_bkgd)

df_level<-data.frame('year'=c(rep('Ctrl_CryEx', nrow(input_CryEx)), rep('KD_CryEx', nrow(input_CryEx)),rep('Ctrl_IR',nrow(input_IR)), rep('KD_IR',nrow(input_IR)),
                              rep('Ctrl_bkgd',nrow(input_bkgd)), rep('KD_bkgd',nrow(input_bkgd))),
                     'MeanRatio'=c(input_CryEx$V5, input_CryEx$V7, input_IR$V5, input_IR$V7, input_bkgd$V5, input_bkgd$V7),
                     'paired'=c(c(1:nrow(input_CryEx)), c(1:nrow(input_CryEx)),
                                c((nrow(input_CryEx)+1):(nrow(input_CryEx)+nrow(input_IR))), c((nrow(input_CryEx)+1):(nrow(input_CryEx)+nrow(input_IR))),
                                c((nrow(input_CryEx)+nrow(input_IR)+1):(nrow(input_CryEx)+nrow(input_IR)+nrow(input_bkgd))), c((nrow(input_CryEx)+nrow(input_IR)+1):(nrow(input_CryEx)+nrow(input_IR)+nrow(input_bkgd))))
                     )

df_level %>%
  ggplot(aes(year,MeanRatio, fill=year))+
  geom_boxplot()+
  ylim(0,1)+
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=year,group=paired), size=0.3,shape=21,position = position_dodge(0.2))+
  theme(legend.position = "none")
wilcox.test(input_bkgd$V7, input_IR$V7)
#############################################
df_bkgd<-data.frame('year'=c(rep('Ctrl_bkgd', nrow(input_bkgd)), rep('KD_bkgd', nrow(input_bkgd))),
                     'MeanRatio'=c(input_bkgd$V5, input_bkgd$V7),
                     'pabkgded'=c(c(1:nrow(input_bkgd)), c(1:nrow(input_bkgd)))
)

df_bkgd %>%
  ggplot(aes(year,MeanRatio, fill=year))+
  geom_boxplot()+
  ylim(0,1)+
  geom_line(aes(group=pabkgded), position = position_dodge(0.2)) +
  geom_point(aes(fill=year,group=pabkgded), size=0.3,shape=21,position = position_dodge(0.2))+
  theme(legend.position = "none")
wilcox.test(input_bkgd$V7, input_bkgd$V7)
################################################################################

numerator<-c((182/612838), (408/6484336), (6720/187149245))
denominator<-c((408+6720)/(6484336+187149245),
               (182+6720)/(612838+187149245),
               (182+408)/(612838+6484336+91128))
barplot(log2(numerator/denominator))
#######################################################################################
numerator<-c((53/872), (52/597),(2148/566847))

denominator<-c((52+6+2148)/(597+712+566847),
               (53+6+2148)/(872+712+566847),
               (53+52+2148)/(872+597+566847),
               (53+52+6)/(872+597+712))
barplot(log2(numerator/denominator))

# RBM34_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235305599 & add_ECL_EKD$Start<235305892,]$MeanRatio.y
# RBM34_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235305961 & add_ECL_EKD$Start<235306262,]$MeanRatio.y
# RBM34_3_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235306769 & add_ECL_EKD$Start<235307050,]$MeanRatio.y
# RBM34_4_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235307215 & add_ECL_EKD$Start<235307487,]$MeanRatio.y
# RBM34_5_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235307919 & add_ECL_EKD$Start<235308072,]$MeanRatio.y
# RBM34_6_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235308076 & add_ECL_EKD$Start<235308385,]$MeanRatio.y
# RBM34_7_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235308597 & add_ECL_EKD$Start<235308892,]$MeanRatio.y
# RBM34_8_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235309105 & add_ECL_EKD$Start<235309267,]$MeanRatio.y
# RBM34_9_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235309933 & add_ECL_EKD$Start<235310235,]$MeanRatio.y
# RBM34_10_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235309282 & add_ECL_EKD$Start<235309560,]$MeanRatio.y
# RBM34_11_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>235310272 & add_ECL_EKD$Start<235310580,]$MeanRatio.y
# 

boxplot(c(RBM34_1, RBM34_2, RBM34_3, RBM34_4, RBM34_5, RBM34_6, RBM34_7, RBM34_8, RBM34_9, RBM34_10, RBM34_11), 
        c(RBM34_1_KD, RBM34_2_KD, RBM34_3_KD, RBM34_4_KD, RBM34_5_KD, RBM34_6_KD, RBM34_7_KD, RBM34_8_KD, RBM34_9_KD, RBM34_10_KD, RBM34_11_KD))
wilcox.test(c(RBM34_1, RBM34_2, RBM34_3, RBM34_4, RBM34_5, RBM34_6, RBM34_7, RBM34_8, RBM34_9, RBM34_10, RBM34_11), 
            c(RBM34_1_KD, RBM34_2_KD, RBM34_3_KD, RBM34_4_KD, RBM34_5_KD, RBM34_6_KD, RBM34_7_KD, RBM34_8_KD, RBM34_9_KD, RBM34_10_KD, RBM34_11_KD))$p.value
MED15_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr22" & add_ECL_EKD$Start>20927249 & add_ECL_EKD$Start<20927554,]$MeanRatio.x
MED15_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr22" & add_ECL_EKD$Start>20926927 & add_ECL_EKD$Start<20927229,]$MeanRatio.x

MED15_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr22" & add_ECL_EKD$Start>20927249 & add_ECL_EKD$Start<20927554,]$MeanRatio.y
MED15_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr22" & add_ECL_EKD$Start>20926927 & add_ECL_EKD$Start<20927229,]$MeanRatio.y

boxplot(c(MED15_1, MED15_2), c(MED15_1_KD, MED15_2_KD))

MAN1B1_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr9" & add_ECL_EKD$Start>139997675 & add_ECL_EKD$Start<139998205,]$MeanRatio.x
MAN1B1_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr9" & add_ECL_EKD$Start>139996230 & add_ECL_EKD$Start<139996725,]$MeanRatio.x

MAN1B1_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr9" & add_ECL_EKD$Start>139997675 & add_ECL_EKD$Start<139998205,]$MeanRatio.y
MAN1B1_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr9" & add_ECL_EKD$Start>139996230 & add_ECL_EKD$Start<139996725,]$MeanRatio.y

boxplot(c(MAN1B1_1, MAN1B1_2), c(MAN1B1_1_KD, MAN1B1_2_KD))

RPRD1B_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr20" & add_ECL_EKD$Start>36706296 & add_ECL_EKD$Start<36706584,]$MeanRatio.x
RPRD1B_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr20" & add_ECL_EKD$Start>36710785 & add_ECL_EKD$Start<36711088,]$MeanRatio.x

RPRD1B_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr20" & add_ECL_EKD$Start>36706296 & add_ECL_EKD$Start<36706584,]$MeanRatio.y
RPRD1B_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr20" & add_ECL_EKD$Start>36710785 & add_ECL_EKD$Start<36711088,]$MeanRatio.y

boxplot(c(RPRD1B_1, RPRD1B_2), c(RPRD1B_1_KD, RPRD1B_2_KD))

TRAPPC10_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr21" & add_ECL_EKD$Start>45461255 & add_ECL_EKD$Start<45461551,]$MeanRatio.x
TRAPPC10_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr21" & add_ECL_EKD$Start>45459466 & add_ECL_EKD$Start<45459761,]$MeanRatio.x

TRAPPC10_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr21" & add_ECL_EKD$Start>45461255 & add_ECL_EKD$Start<45461551,]$MeanRatio.y
TRAPPC10_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr21" & add_ECL_EKD$Start>45459466 & add_ECL_EKD$Start<45459761,]$MeanRatio.y

boxplot(c(TRAPPC10_1, TRAPPC10_2), c(TRAPPC10_1_KD, TRAPPC10_2_KD))
wilcox.test(c(TRAPPC10_1, TRAPPC10_2), c(TRAPPC10_1_KD, TRAPPC10_2_KD))

LRP11_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr6" & add_ECL_EKD$Start>150178386 & add_ECL_EKD$Start<150178670,]$MeanRatio.x
LRP11_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr6" & add_ECL_EKD$Start>150180534 & add_ECL_EKD$Start<150180831,]$MeanRatio.x

LRP11_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr6" & add_ECL_EKD$Start>150178386 & add_ECL_EKD$Start<150178670,]$MeanRatio.y
LRP11_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr6" & add_ECL_EKD$Start>150180534 & add_ECL_EKD$Start<150180831,]$MeanRatio.y

boxplot(c(LRP11_1, LRP11_2), c(LRP11_1_KD, LRP11_2_KD))
wilcox.test(c(LRP11_1, LRP11_2), c(LRP11_1_KD, LRP11_2_KD))

AUH_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr9" & add_ECL_EKD$Start>94100714 & add_ECL_EKD$Start<94100984,]$MeanRatio.x
AUH_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr9" & add_ECL_EKD$Start>94098632 & add_ECL_EKD$Start<94098933,]$MeanRatio.x

AUH_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr9" & add_ECL_EKD$Start>94100714 & add_ECL_EKD$Start<94100984,]$MeanRatio.y
AUH_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr9" & add_ECL_EKD$Start>94098632 & add_ECL_EKD$Start<94098933,]$MeanRatio.y

boxplot(c(AUH_1, AUH_2), c(AUH_1_KD, AUH_2_KD))
wilcox.test(c(AUH_1, AUH_2), c(AUH_1_KD, AUH_2_KD))

ZFP62_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr5" & add_ECL_EKD$Start>180283673 & add_ECL_EKD$Start<180283978,]$MeanRatio.x
ZFP62_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr5" & add_ECL_EKD$Start>180284491 & add_ECL_EKD$Start<180284794,]$MeanRatio.x

ZFP62_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr5" & add_ECL_EKD$Start>180283673 & add_ECL_EKD$Start<180283978,]$MeanRatio.y
ZFP62_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr5" & add_ECL_EKD$Start>180284491 & add_ECL_EKD$Start<180284794,]$MeanRatio.y

boxplot(c(ZFP62_1, ZFP62_2), c(ZFP62_1_KD, ZFP62_2_KD))
wilcox.test(c(ZFP62_1, ZFP62_2), c(ZFP62_1_KD, ZFP62_2_KD))


RAP1GAP2_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr17" & add_ECL_EKD$Start>2924723 & add_ECL_EKD$Start<2925028,]$MeanRatio.x
RAP1GAP2_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr17" & add_ECL_EKD$Start>2928389 & add_ECL_EKD$Start<2928677,]$MeanRatio.x

RAP1GAP2_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr17" & add_ECL_EKD$Start>2924723 & add_ECL_EKD$Start<2925028,]$MeanRatio.y
RAP1GAP2_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr17" & add_ECL_EKD$Start>2928389 & add_ECL_EKD$Start<2928677,]$MeanRatio.y

boxplot(c(RAP1GAP2_1, RAP1GAP2_2), c(RAP1GAP2_1_KD, RAP1GAP2_1_KD))
wilcox.test(c(RAP1GAP2_1, RAP1GAP2_2), c(RAP1GAP2_1_KD, RAP1GAP2_2_KD))

SLC2A9_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr4" & add_ECL_EKD$Start>9871458 & add_ECL_EKD$Start<9871881,]$MeanRatio.x
SLC2A9_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr4" & add_ECL_EKD$Start>9872913 & add_ECL_EKD$Start<9873363,]$MeanRatio.x

SLC2A9_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr4" & add_ECL_EKD$Start>9871458 & add_ECL_EKD$Start<9871881,]$MeanRatio.y
SLC2A9_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr4" & add_ECL_EKD$Start>9872913 & add_ECL_EKD$Start<9873363,]$MeanRatio.y

boxplot(c(SLC2A9_1, SLC2A9_2), c(SLC2A9_1_KD, SLC2A9_1_KD))
wilcox.test(c(SLC2A9_1, SLC2A9_2), c(SLC2A9_1_KD, SLC2A9_2_KD))



DENND6A_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr3" & add_ECL_EKD$Start>57636513 & add_ECL_EKD$Start<57636810,]$MeanRatio.x
DENND6A_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr3" & add_ECL_EKD$Start>57637654 & add_ECL_EKD$Start<57637774,]$MeanRatio.x
DENND6A_3<-add_ECL_EKD[add_ECL_EKD$chr=="chr3" & add_ECL_EKD$Start>57639191 & add_ECL_EKD$Start<57639490,]$MeanRatio.x
DENND6A_4<-add_ECL_EKD[add_ECL_EKD$chr=="chr3" & add_ECL_EKD$Start>57639669 & add_ECL_EKD$Start<57639981,]$MeanRatio.x

DENND6A_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr3" & add_ECL_EKD$Start>57636513 & add_ECL_EKD$Start<57636810,]$MeanRatio.y
DENND6A_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr3" & add_ECL_EKD$Start>57637654 & add_ECL_EKD$Start<57637774,]$MeanRatio.y
DENND6A_3_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr3" & add_ECL_EKD$Start>57639191 & add_ECL_EKD$Start<57639490,]$MeanRatio.y
DENND6A_4_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr3" & add_ECL_EKD$Start>57639669 & add_ECL_EKD$Start<57639981,]$MeanRatio.y

boxplot(c(DENND6A_1, DENND6A_2, DENND6A_3, DENND6A_4), c(DENND6A_1_KD, DENND6A_2_KD, DENND6A_3_KD, DENND6A_4_KD))
wilcox.test(c(DENND6A_1, DENND6A_2, DENND6A_3, DENND6A_4), c(DENND6A_1_KD, DENND6A_2_KD, DENND6A_3_KD, DENND6A_4_KD))

HP1BP3_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>21081779 & add_ECL_EKD$Start<21082081,]$MeanRatio.x
HP1BP3_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>21083202 & add_ECL_EKD$Start<2108350,]$MeanRatio.x
HP1BP3_3<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>21082929 & add_ECL_EKD$Start<21083070,]$MeanRatio.x

HP1BP3_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>21081779 & add_ECL_EKD$Start<21082081,]$MeanRatio.y
HP1BP3_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>21083202 & add_ECL_EKD$Start<2108350,]$MeanRatio.y
HP1BP3_3_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>21082929 & add_ECL_EKD$Start<21083070,]$MeanRatio.y

boxplot(c(HP1BP3_1, HP1BP3_2, HP1BP3_3), c(HP1BP3_1_KD, HP1BP3_2_KD, HP1BP3_2_KD))
wilcox.test(c(HP1BP3_1, HP1BP3_2, HP1BP3_3), c(HP1BP3_1_KD, HP1BP3_2_KD, HP1BP3_2_KD))

NUDT19_1<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33191429 & add_ECL_EKD$Start<33191502,]$MeanRatio.x
NUDT19_2<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33191503 & add_ECL_EKD$Start<33191636,]$MeanRatio.x
NUDT19_3<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33191669 & add_ECL_EKD$Start<33191960,]$MeanRatio.x
NUDT19_4<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33194557 & add_ECL_EKD$Start<33194868,]$MeanRatio.x
NUDT19_5<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33195172 & add_ECL_EKD$Start<33195465,]$MeanRatio.x
NUDT19_6<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33195731 & add_ECL_EKD$Start<33196034,]$MeanRatio.x
NUDT19_7<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33192106 & add_ECL_EKD$Start<33192403,]$MeanRatio.x
NUDT19_8<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33192776 & add_ECL_EKD$Start<33192907,]$MeanRatio.x
NUDT19_9<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33193410 & add_ECL_EKD$Start<33193530,]$MeanRatio.x
NUDT19_10<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33193914 & add_ECL_EKD$Start<33194206,]$MeanRatio.x
NUDT19_11<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33194262 & add_ECL_EKD$Start<33194554,]$MeanRatio.x

NUDT19_1_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33191429 & add_ECL_EKD$Start<33191502,]$MeanRatio.y
NUDT19_2_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33191503 & add_ECL_EKD$Start<33191636,]$MeanRatio.y
NUDT19_3_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33191669 & add_ECL_EKD$Start<33191960,]$MeanRatio.y
NUDT19_4_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33194557 & add_ECL_EKD$Start<33194868,]$MeanRatio.y
NUDT19_5_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33195172 & add_ECL_EKD$Start<33195465,]$MeanRatio.y
NUDT19_6_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33195731 & add_ECL_EKD$Start<33196034,]$MeanRatio.y
NUDT19_7_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33192106 & add_ECL_EKD$Start<33192403,]$MeanRatio.y
NUDT19_8_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33192776 & add_ECL_EKD$Start<33192907,]$MeanRatio.y
NUDT19_9_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33193410 & add_ECL_EKD$Start<33193530,]$MeanRatio.y
NUDT19_10_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33193914 & add_ECL_EKD$Start<33194206,]$MeanRatio.y
NUDT19_11_KD<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33194262 & add_ECL_EKD$Start<33194554,]$MeanRatio.y

boxplot(c(NUDT19_1, NUDT19_2, NUDT19_3, NUDT19_4, NUDT19_5, NUDT19_6, NUDT19_7, NUDT19_8, NUDT19_9, NUDT19_10, NUDT19_11), 
        c(NUDT19_1_KD, NUDT19_2_KD, NUDT19_3_KD, NUDT19_4_KD, NUDT19_5_KD, NUDT19_6_KD, NUDT19_7_KD, NUDT19_8_KD, NUDT19_9_KD, NUDT19_10_KD, NUDT19_11_KD))
wilcox.test(c(NUDT19_1, NUDT19_2, NUDT19_3, NUDT19_4, NUDT19_5, NUDT19_6, NUDT19_7, NUDT19_8, NUDT19_9, NUDT19_10, NUDT19_11), 
            c(NUDT19_1_KD, NUDT19_2_KD, NUDT19_3_KD, NUDT19_4_KD, NUDT19_5_KD, NUDT19_6_KD, NUDT19_7_KD, NUDT19_8_KD, NUDT19_9_KD, NUDT19_10_KD, NUDT19_11_KD))$p.value



TRAPPC10<-add_ECL_EKD[add_ECL_EKD$chr=="chr21" & add_ECL_EKD$Start>45458057 & add_ECL_EKD$Start<45462429,]$Site_ID
SLC2A9<-add_ECL_EKD[add_ECL_EKD$chr=="chr4" & add_ECL_EKD$Start>9870453 & add_ECL_EKD$Start<9876758,]$Site_ID
ZFP62<-add_ECL_EKD[add_ECL_EKD$chr=="chr5" & add_ECL_EKD$Start>180283138 & add_ECL_EKD$Start<180285053,]$Site_ID
MAN1B1<-add_ECL_EKD[add_ECL_EKD$chr=="chr9" & add_ECL_EKD$Start>139996989 & add_ECL_EKD$Start<140000368,]$Site_ID
LRP11<-add_ECL_EKD[add_ECL_EKD$chr=="chr6" & add_ECL_EKD$Start>150178212 & add_ECL_EKD$Start<150139893,]$Site_ID
DENND6A<-add_ECL_EKD[add_ECL_EKD$chr=="chr3" & add_ECL_EKD$Start>57636367 & add_ECL_EKD$Start<57639981,]$Site_ID
NUDT19<-add_ECL_EKD[add_ECL_EKD$chr=="chr19" & add_ECL_EKD$Start>33183580 & add_ECL_EKD$Start<33200090,]$Site_ID
RAP1GAP2<-add_ECL_EKD[add_ECL_EKD$chr=="chr17" & add_ECL_EKD$Start>2924723 & add_ECL_EKD$Start<2928677,]$Site_ID
RPRD1B<-add_ECL_EKD[add_ECL_EKD$chr=="chr20" & add_ECL_EKD$Start>36706296 & add_ECL_EKD$Start<36709226,]$Site_ID
HP1BP3<-add_ECL_EKD[add_ECL_EKD$chr=="chr1" & add_ECL_EKD$Start>21081779 & add_ECL_EKD$Start<21083500,]$Site_ID
AUH<-add_ECL_EKD[add_ECL_EKD$chr=="chr9" & add_ECL_EKD$Start>94098345 & add_ECL_EKD$Start<94105657,]$Site_ID

###### Highlight dots ######
library(dplyr)
library(ggplot2)
ECL_EKD %>%
  ggplot(aes(x=MeanRatio.x, y=MeanRatio.y))+
  geom_point(alpha=0.3)
tot<-c(TRAPPC10, SLC2A9, ZFP62, MAN1B1, LRP11, DENND6A, NUDT19, RAP1GAP2, HP1BP3, AUH)
highlight_df<-data.frame()
for (k in 1:length(tot)){
  highlight_df<-rbind(highlight_df, ECL_EKD[ECL_EKD$Site_ID==tot[k],])
}

ECL_EKD %>%
  ggplot(aes(x=MeanRatio.x, y=MeanRatio.y))+
  geom_point(alpha=0.3)+
  geom_point(data = highlight_df,
             aes(x=MeanRatio.x, y=MeanRatio.y),
             color='red',
             size=3)

################################################################################
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
# HNRNPM_RAW_Counts/
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
#write.table(hnM_tot, '~/Downloads/hnM_tot.txt', sep = '\t', quote=F, row.names = F)
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),
  values=hnM_tot$gene,
  mart=mart)
head(genes)
dim(genes)
symbol_HMLE<-merge(genes, hnM_tot, by.x="ensembl_gene_id", by.y="gene")
symbol_HMLE<-symbol_HMLE[!duplicated(symbol_HMLE$ensembl_gene_id),]
head(symbol_HMLE)
#write.table(symbol_HMLE, '~/Downloads/hnM_tot.txt', sep = '\t', quote=F, row.names = F)
#hnM_tot_tmp<-read.table('~/Downloads/hnM_tot.txt', sep = '\t', header = T)
#hnM_tot<-hnM_tot_tmp[,-c(2:3)]
hnM_tot<-symbol_HMLE[,-c(2:3)]
colnames(hnM_tot)<-c('gene','len','HTE1','HTE2','HTE3','HTE4','HH5','HH6')
#colnames(hnM_tot)<-c('gene','len','HTE1','HTE2','HTE3','HTE4','HH_5','HH_6','HH5','HH6')
#### get rid of low abundance genes
#sum_lib<-c(sum(hnM_tot$HTE1), sum(hnM_tot$HTE2), sum(hnM_tot$HTE3), sum(hnM_tot$HTE4), sum(hnM_tot$HH_5), sum(hnM_tot$HH_6), sum(hnM_tot$HH5), sum(hnM_tot$HH6))
sum_lib<-c(sum(hnM_tot$HTE1), sum(hnM_tot$HTE2), sum(hnM_tot$HTE3), sum(hnM_tot$HTE4), sum(hnM_tot$HH5), sum(hnM_tot$HH6))
keep <- rowSums(cpm(hnM_tot[,-c(1:2)])>(5*1000000/(min(sum_lib))))>=3
#hnM_keep <- hnM_tot[keep,]
head(hnM_keep)
rownames(hnM_keep)<-hnM_keep$gene
# 29221
dim(hnM_keep)



corrected_data<-data.frame(corrected_data)
head(corrected_data)

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
write.table(symbol_HMLE$hgnc_symbol, '~/Downloads/input_background_file.txt', sep = '\t', row.names = F, quote = F)
### rpkm_file: 
# count up the total reads in a sample and divide that number by 1000000
head(symbol_HMLE)
symbol_HMLE<-symbol_HMLE[!duplicated(symbol_HMLE$hgnc_symbol),]
dim(symbol_HMLE)

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
dim(convert_HTE)
final_df<-convert_HTE[,c(2:7)]
head(final_df)
### 
for (l in 1:nrow(final_df)){
  final_df[l,]<-(as.numeric(convert_HTE[l,2:7])*1000)/convert_HTE[l,9]
}
rownames(final_df)<-convert_HTE$ensembl
head(final_df)
write.table(data.frame('gene'=rownames(final_df), 'RPKM'=final_df$HTE1), '~/Downloads/rpkm_file.txt', sep = '\t', row.names = F, quote = F)

################################################################################
summary(final_df$HTE1)

input_keep<-convert_HTE[,c(2:7)]
rownames(input_keep)<-convert_HTE$V1
head(input_keep)
summary(input_keep)




################################################################################
kp_y<-DGEList(cor_df)
kp_y$samples
kp_y<-calcNormFactors(kp_y, method='TMM')
lcpmtmm_hnM<-cpm(kp_y, log=F, normalized.lib.sizes=T)
dim(lcpmtmm_hnM)
head(lcpmtmm_hnM)
# convert Ensembl gene id to entrez ID
#library(biomaRt)
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),
  values=rownames(lcpmtmm_hnM),
  mart=mart)
head(genes)
dim(genes)
symbol_HMLE<-merge(genes, data.frame(cbind('ensembl_gene_id'=rownames(lcpmtmm_hnM), lcpmtmm_hnM)), by="ensembl_gene_id")
#symbol_HMLE<-symbol_HMLE[!duplicated(symbol_HMLE$ensembl_gene_id),]
symbol_HMLE<-symbol_HMLE[!duplicated(symbol_HMLE$hgnc_symbol),]
head(symbol_HMLE)
rownames(symbol_HMLE)<-symbol_HMLE$hgnc_symbol
#write.table(symbol_HMLE, '~/Downloads/symbo_hmle.txt', row.names = F, sep = '\t', quote = F)
#symbol_HMLE<-read.table('~/Downloads/symbo_hmle.txt', header = T, sep = '\t')
#########batch correct############
# conditions = pheno$cancer
# library_methods = pheno$batch
# groups = sapply(as.character(conditions), switch, "Ctrl"=1, "KD"=2, USE.NAMES = F)
# batches = as.numeric(pheno$batch)
# corrected_data = ComBat_seq(counts = as.matrix(hnM_keep[,-c(1:2)]),
#                             batch = batches, group = groups)
# pca_corrected_obj = prcomp(hnM_keep[,-c(1:2)])
# pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)
# pca_corrected[,"condition"] = conditions
# pca_corrected[,"library_method"] = as.factor(library_methods)
# 
# cols<-c("Ctrl"="#481567FF", "KD"="#1F968BFF")
# p2 = ggplot(data=pca_corrected, aes(x=PC1, y=PC2, color=condition, shape=library_method))
# p2 = p2 + geom_point(size=3)
# p2 = p2 + stat_ellipse(type = "norm", linetype=2)
# p2 = p2 + labs(title="PCA, RNA-seq counts for 16 HBR/UHR and Ribo/PolyA samples (batch corrected data)", color="Condition", shape="Library Method")
# p2 = p2 + scale_colour_manual(values = cols)
# 
# pdf(file="UncorrectedPCA.pdf")
# library(gridExtra)
# library(ggplot2)
# grid.arrange(p1, p2, nrow = 2)
# dev.off()
# ##################################
# BiocManager::install('bladderbatch')
# library(bladderbatch)
# data(bladderdata)
# 
# dat<-bladderEset[1:50,]
# pheno = pData(dat)
# edata = exprs(dat)
# batch = pheno$batch
# mod = model.matrix(~as.factor(cancer), data=pheno)
# 
# combat_edata2 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, prior.plots = FALSE)
# BiocManager::install('sva')
# library(sva)
# library(edgeR)
# library(UpSetR)
# batch = c(1,1,1,1,2,2,3,3)
# pheno<-data.frame('sample'=c(1,2,3,4,5,6,7,8), 'batch'=batch, 'cancer'=c('Ctrl','KD','Ctrl','KD','Ctrl','KD','Ctrl','KD'))
# rownames(pheno)<-colnames(symbol_HMLE[,-c(1:3)])
# symbol_HMLE<-na.omit(symbol_HMLE[!duplicated(symbol_HMLE$hgnc_symbol),])
# rownames(symbol_HMLE)<-symbol_HMLE$hgnc_symbol
# for (l in 4:ncol(symbol_HMLE)){
#   symbol_HMLE[,l]<-as.numeric(symbol_HMLE[,l])
# }
# combat_data<-ComBat(dat=symbol_HMLE[,-c(1:3)], batch=batch, mod=NULL, par.prior=FALSE, mean.only = TRUE)
# combat_data_new<-ComBat(dat=symbol_HMLE[,-c(1:3)], batch=batch, mod=mod, par.prior=TRUE, prior.plots = FALSE)
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
head(hnM_keep[,-c(1:2)])
# pca_corrected_obj = prcomp(hnM_keep[,-c(1:2)])
# pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)
# pca_corrected[,"condition"] = conditions
# pca_corrected[,"library_method"] = as.factor(library_methods)
# 
# cols<-c("Ctrl"="#481567FF", "KD"="#1F968BFF")
# p2 = ggplot(data=pca_corrected, aes(x=PC1, y=PC2, color=condition, shape=library_method))
# p2 = p2 + geom_point(size=3)
# p2 = p2 + stat_ellipse(type = "norm", linetype=2)
# p2 = p2 + labs(title="PCA, RNA-seq counts for 16 HBR/UHR and Ribo/PolyA samples (batch corrected data)", color="Condition", shape="Library Method")
# p2 = p2 + scale_colour_manual(values = cols)
# 
#  cat hallmark_interferon_alpha_response_gene_set.txt hallmark_interferon_gamma_response_gene_set.txt A.J.Scadden_nature_struct_bio.interferon_I_list|awk '!seen[$0]++' > ISG.txt
ISG<-read.table('~/Downloads/ISG.txt')
head(ISG)
combat_df<-data.frame('hgnc_symbol'=symbol_HMLE$hgnc_symbol, symbol_HMLE[,-c(1:3)])
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
  FC_input[k,1]<-(ISG_input[k,1]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  FC_input[k,2]<-(ISG_input[k,2]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  FC_input[k,3]<-(ISG_input[k,3]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  FC_input[k,4]<-(ISG_input[k,4]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  FC_input[k,5]<-(ISG_input[k,5]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  FC_input[k,6]<-(ISG_input[k,6]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5)]))+0.01)
  #FC_input[k,7]<-(ISG_input[k,7]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5,7)]))+0.01)
  #FC_input[k,8]<-(ISG_input[k,8]+0.01)/(mean(as.numeric(ISG_input[k,c(1,3,5,7)]))+0.01)
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
### ??? What is the method they use in Fig6C ??? 
library(pheatmap)
library(RColorBrewer)
#install.packages('viridis')
library(viridis)
genelist<-c("MX1","IFI27","IFIT1","IFIT2","IFI3","DDX58","IFIH1","TLR3","CXCL10","IFI6",
            "OAS3","IRF7","ISG15","BST2")
labels<-rownames(rev_rev_trim)
labels[!labels %in% genelist] <- ""
#breaksList = seq(0, 100, by = 1)
# [1:75,]

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
################################################################################
###########################GSEA########################
BiocManager::install('fgsea')
library(fgsea)
library(ggplot2)
# cat c2.all.v7.4.symbols.gmt.txt|grep "REACTOME"> c2_REACTOME.gmt.txt
#filename = "~/Downloads/c2_REACTOME.gmt.txt"
filename = "~/Downloads/c2_REACTOME.gmt.txt"
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
  
  one_pathway = one_pathway[-1]
  
  all_pathway[[i]] = one_pathway
  
}

names(all_pathway) = names
#######################
##################################################
#BiocManager::install('DESeq2')
library('DESeq2')
#BiocManager::install('edgeR')
library('edgeR')


# load in the conversion table:
convert<-read.table('~/Downloads/genesym_to_ensembl_genelengths.txt', header = F, sep = '\t')
head(convert)
# dimensions: 45471;3
dim(convert)
convert_HTE<- merge(data.frame('ensembl'=rownames(cor_df), cor_df), convert, by.x='ensembl', by.y='V2')
head(convert_HTE)
convert_HTE<-convert_HTE[!duplicated(convert_HTE$V1),]
# dimensions: 45471;8
dim(convert_HTE)
input_keep<-convert_HTE[,c(2:7)]
rownames(input_keep)<-convert_HTE$V1
head(input_keep)
summary(input_keep)
########################################
########################################

conditions<-factor(c('Ctrl','KD','Ctrl','KD','Ctrl','KD'))
colData<-data.frame(sampleNames=c('HTE1','HTE2','HTE3','HTE4','HH5','HH6'), conditions=conditions)
row.names(colData)<-colData$sampleNames
dds<-DESeqDataSetFromMatrix(countData = as.matrix(input_keep),
                            colData = colData,
                            design= ~conditions)
dds<-DESeq(dds)
resultsNames(dds)
res<-results(dds, name="conditions_KD_vs_Ctrl")
summary(res)
grp.mean<-sapply(levels(dds$conditions), function(lvl) rowMeans(counts(dds, normalized=TRUE)[,dds$conditions==lvl]))
norm.counts<-counts(dds, normalized=TRUE)
all<-data.frame(res, grp.mean, norm.counts)
head(all)
dim(all)
all_nonNA<-data.frame(all[!is.na(all$log2FoldChange),])
# 27136
dim(all_nonNA)
head(all_nonNA)
#########################
UP_DEG<-all_nonNA[all_nonNA$log2FoldChange> 1.5 & all_nonNA$padj<0.05, c(2,5:6)]
dim(UP_DEG)
head(UP_DEG)

write.table(rownames(UP_DEG), '~/Downloads/UP_DEG.txt', sep = '\t', row.names = F, quote = F)
DN_DEG<-all_nonNA[all_nonNA$log2FoldChange< -1.5 & all_nonNA$padj<0.05, c(2,5:6)]
dim(DN_DEG)
head(DN_DEG)

all_input<-data.frame('gene'=rownames(all_nonNA), 'lfc'=all_nonNA$log2FoldChange)
# Getting the protein coding genes only:
library(biomaRt)
mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
genes<-biomaRt::getBM(attributes = c("external_gene_name","chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"), values = list("protein_coding", c(1:22,"X","Y")), mart = mart)
PC_input<-merge(all_input, genes, by.x="gene", by.y="external_gene_name")
PC_input<-PC_input[!duplicated(PC_input$gene),]
dim(PC_input)
head(PC_input)
all_input<-PC_input[,c(1:2)]

UP_DEG<-all_input[all_input$lfc>1, ]
DN_DEG<-all_input[all_input$lfc< -1, ]
write.table(UP_DEG$gene, '~/Downloads/UP_DEG.txt', sep = '\t', row.names = F, quote = F)

write.table(all_input, '~/Downloads/all_input.txt', sep = '\t', row.names = F, quote = F)

write.table(rownames(DN_DEG), '~/Downloads/DN_DEG.txt', sep = '\t', row.names = F, quote = F)


#########################
input_HMLE<-as.numeric(all_nonNA$log2FoldChange)
names(input_HMLE)<-rownames(all_nonNA)
summary(input_HMLE)
barplot(sort(input_HMLE, decreasing=T))
fgseaRes_HMLE <- fgsea(pathways=all_pathway, 
                       stats=input_HMLE,
                       minSize=15,
                       maxSize=500,
                       nperm=10000)

dim(fgseaRes_HMLE[fgseaRes_HMLE$padj<0.05,])
head(fgseaRes_HMLE[order(padj, -abs(NES)), ], n=30)
neg<-fgseaRes_HMLE[fgseaRes_HMLE$NES<0,]
dim(neg)
head(neg[order(padj, -abs(NES)), ], n=30)
# SEP14
write.table(fgseaRes_HMLE[fgseaRes_HMLE$padj<0.05,c(1,3,5,6,7)],'~/Downloads/HMLE_REACTOME_rep3_OCT22_2022.txt', sep='\t', row.names = FALSE)
##################
input_bar<-read.table('~/Downloads/HMLE_REACTOME_rep3_SEP14_2022.txt', sep = '\t', header = T)
head(input_bar)
dim(input_bar)
order_input_bar<-input_bar[rev(order(input_bar$NES)),]
head(order_input_bar)
barplot(order_input_bar$NES, horiz = T)


plotEnrichment(all_pathway[["REACTOME_ANTIMICROBIAL_PEPTIDES"]],
               input_HMLE)+labs(title="")
#topPathwaysUp<-fgseaRes_HTE[ES>0][head(order(pval),n=10), pathway]
top<-c("REACTOME_INTERFERON_ALPHA_BETA_SIGNALING")
plotGseaTable(all_pathway[top],input_HMLE,
              fgseaRes_HMLE, gseaParam = 0.5)
################################################################################
################################################################################
################################################################################
# Clusterprofile 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

#organism="org.Hs.eg.db"
#BiocManager::install(organism, character.only=TRUE)
#library(organism, character.only=TRUE)

#deg_all<-all_nonNA[all_nonNA$log2FoldChange> 4 & all_nonNA$padj< 0.01, ]
#deg_all<-na.omit(deg_all)
#dim(deg_all)
#head(deg_all)
#library(biomaRt)
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="hgnc_symbol",
  attributes=c("hgnc_symbol","entrezgene_id","ensembl_gene_id"),
  values=rownames(all_nonNA),
  mart=mart)
head(genes)
dim(genes)
symbol_brca<-merge(genes, data.frame(cbind("hgnc_symbol"=rownames(all_nonNA),all_nonNA)), by="hgnc_symbol")
symbol_brca<-symbol_brca[!duplicated(symbol_brca$hgnc_symbol),]
head(symbol_brca)
# 1449
dim(symbol_brca)
deg_all<-symbol_brca[symbol_brca$log2FoldChange>2 & symbol_brca$padj<0.05,]
head(deg_all)
dim(deg_all)

# original_gene_list<-as.numeric(symbol_brca$log2FoldChange)
# names(original_gene_list)<-symbol_brca$entrezgene_id
# gene_list<-na.omit(original_gene_list)
# gene_list=sort(gene_list, decreasing = TRUE)
# 
### enrich KEGG###
KEGG_enrich<-enrichKEGG(deg_all$entrezgene_id, organism="hsa", keyType="kegg", pvalueCutoff=0.05, pAdjustMethod="BH",
                        universe=as.character(symbol_brca$entrezgene_id), minGSSize=5, maxGSSize=500, qvalueCutoff=0.2, use_internal_data=FALSE)

summary(KEGG_enrich)
### GO ###

### GSEA ###

gse<-enrichGO(gene = deg_all$entrezgene_id,
              ont = "BP",
              universe= as.character(symbol_brca$entrezgene_id),
              #keyType = "ENTREZID",
              #nPerm=10000,
              # minGSSize = 5,
              #maxGSSize = 800,
              #pvalueCutoff = 0.05,
              pAdjustMethod="BH",
              #verbose = TRUE,
              OrgDb = org.Hs.eg.db,
              qvalueCutoff = 0.01)
#require(DOSE)
dotplot(gse, showCategory=100, split=".sign")+facet_grid(.~.sign)
emapplot(gse, showCategory = 10)




















pca_corrected_obj = prcomp(new_rev_trim)
pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)
pca_corrected[,"condition"] = conditions
pca_corrected[,"library_method"] = as.factor(library_methods)

cols<-c("Ctrl"="#481567FF", "KD"="#1F968BFF")
p2 = ggplot(data=pca_corrected, aes(x=PC1, y=PC2, color=condition, shape=library_method))
p2 = p2 + geom_point(size=3)
p2 = p2 + stat_ellipse(type = "norm", linetype=2)
p2 = p2 + labs(title="PCA, RNA-seq counts for 16 HBR/UHR and Ribo/PolyA samples (batch corrected data)", color="Condition", shape="Library Method")
p2 = p2 + scale_colour_manual(values = cols)

pdf(file="ISG_PCA.pdf")
library(gridExtra)
library(ggplot2)
grid.arrange(p1, p2, nrow = 2)
dev.off()
################################################################################
################################################################################
################################################################################
# cat LINEs.bed|awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}'|tr -d '"'> LINE_trim.bed
# cat SINEs.bed|awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}'|tr -d '"'> SINE_trim.bed
# cat LTR.bed|awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}'|tr -d '"'> LTR_trim.bed
# cat DNA.bed|awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}'|tr -d '"'> DNA_trim.bed
# cat Repeats_all.bed|awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"$9}'|tr -d '"'> Repeats_trim.bed
# awk 'NR==FNR{a[$0];next} !($0 in a)' LINE_trim.bed Repeats_trim.bed > tmp
# awk 'NR==FNR{a[$0];next} !($0 in a)' SINE_trim.bed tmp > tmp1
# awk 'NR==FNR{a[$0];next} !($0 in a)' LTR_trim.bed tmp1 > tmp2
# awk 'NR==FNR{a[$0];next} !($0 in a)' DNA_trim.bed tmp2 > Other_repeats_trim.bed
# cat IR_SIG_hnM.txt|awk '{print $4"\t"$6"\t"$7"\t"$2"\t"$3"\t"$5}'|awk '!seen[$0]++'>IR_SIG_gene.bed
# cat CryEx_M_sig_PSI.txt|awk '{print $4"\t"$6"\t"$7"\t"$2"\t"$3"\t"$5}'|awk '!seen[$0]++'> MIDDLE_gene.bed
# cat annotated_Ctrl_hnM_KD.txt|awk '{print $4"\t"$6"\t"$7"\t"$2"\t"$3"\t"$5}'|awk '!seen[$0]++'> annotated_AS_gene.bed
###################################
# dsRNA
# SINE/Alu & LINE/L1
# Other Alu in it
# >500nt
# select by RNAfold -> narrow down the list 

# cat LINE_trim.bed|awk '{print $1"\t"($2-1)"\t"($3+1)"\t"$4"\t"$5"\t"$6}'> LINE_1_trim.bed 
# cat SINE_trim.bed|awk '{print $1"\t"($2-1)"\t"($3+1)"\t"$4"\t"$5"\t"$6}'> SINE_1_trim.bed
# bedtools intersect -a LINE_1_trim.bed -b SINE_1_trim.bed -wa -wb|awk '!seen[$0]++'|awk '($5=="LINE/L1")'|awk '($11=="SINE/Alu")'> L1_Alu_trim.bed
# bedtools intersect -a IR_SIG_gene.bed -b L1_Alu_trim.bed -wa|awk '($3-$2+1)>500'|awk '!seen[$0]++'> IR_L1_Alu
## bedtools intersect -a LINE_1_trim.bed -b SINE_1_trim.bed -wa -wb|awk '!seen[$0]++'|awk '($5=="LINE/L1")'|awk '($11=="SINE/Alu")'|awk '{print $7"\t"($8+1)"\t"($9-1)"\t"$10"\t"$11"\t"$12}'|awk '!seen[$0]++'> Alu_trim.bed
## grep -v -f Alu_trim.bed SINE_trim.bed > other_SINE_trim.bed

# To see if they are expressed in all included cancers
IR_L1_Alu<-read.table('~/Downloads/MIDDLE_L1_Alu', sep = '\t')
head(IR_L1_Alu)
dim(IR_L1_Alu)
LUSC<-read.table('~/Downloads/lq_hq_0.16_LUSC.txt', sep = '\t')
head(LUSC)
dim(LUSC)
IR_LUSC<-merge(IR_L1_Alu, data.frame('genes'=rownames(LUSC), LUSC), by.x = 'V5', by.y = 'genes')
IR_LUSC<-IR_LUSC[!duplicated(IR_LUSC$V5),c(1,5)]
dim(IR_LUSC)
head(IR_LUSC)
colnames(IR_LUSC)<-c('genes','LUSC')
TGCT<-read.table('~/Downloads/lq_hq_0.16_TGCT.txt', sep = '\t')
IR_TGCT<-merge(IR_L1_Alu, data.frame('genes'=rownames(TGCT), TGCT), by.x = 'V5', by.y = 'genes')
IR_TGCT<-IR_TGCT[!duplicated(IR_TGCT$V5),c(1,5)]
dim(IR_TGCT)
head(IR_TGCT)
colnames(IR_TGCT)<-c('genes','TGCT')
OV<-read.table('~/Downloads/lq_hq_0.16_OV.txt', sep = '\t')
IR_OV<-merge(IR_L1_Alu, data.frame('genes'=rownames(OV), OV), by.x = 'V5', by.y = 'genes')
IR_OV<-IR_OV[!duplicated(IR_OV$V5),c(1,5)]
dim(IR_OV)
head(IR_OV)
colnames(IR_OV)<-c('genes','OV')
GBM<-read.table('~/Downloads/lq_hq_0.16_GBM.txt', sep = '\t')
IR_GBM<-merge(IR_L1_Alu, data.frame('genes'=rownames(GBM), GBM), by.x = 'V5', by.y = 'genes')
IR_GBM<-IR_GBM[!duplicated(IR_GBM$V5),c(1,5)]
dim(IR_GBM)
head(IR_GBM)
colnames(IR_GBM)<-c('genes','GBM')
DLBC<-read.table('~/Downloads/lq_hq_0.16_DLBC.txt', sep = '\t')
IR_DLBC<-merge(IR_L1_Alu, data.frame('genes'=rownames(DLBC), DLBC), by.x = 'V5', by.y = 'genes')
IR_DLBC<-IR_DLBC[!duplicated(IR_DLBC$V5),c(1,5)]
dim(IR_DLBC)
head(IR_DLBC)
colnames(IR_DLBC)<-c('genes','DLBC')
LUAD<-read.table('~/Downloads/lq_hq_0.16_LUAD.txt', sep = '\t')
IR_LUAD<-merge(IR_L1_Alu, data.frame('genes'=rownames(LUAD), LUAD), by.x = 'V5', by.y = 'genes')
IR_LUAD<-IR_LUAD[!duplicated(IR_LUAD$V5),c(1,5)]
dim(IR_LUAD)
head(IR_LUAD)
colnames(IR_LUAD)<-c('genes','LUAD')
BASAL<-read.table('~/Downloads/lq_hq_0.16_Basal.txt', sep = '\t')
IR_BASAL<-merge(IR_L1_Alu, data.frame('genes'=rownames(BASAL), BASAL), by.x = 'V5', by.y = 'genes')
IR_BASAL<-IR_BASAL[!duplicated(IR_BASAL$V5),c(1,5)]
dim(IR_BASAL)
head(IR_BASAL)
colnames(IR_BASAL)<-c('genes','BASAL')
IR_tot<-merge(merge(merge(merge(merge(merge(IR_LUSC, IR_TGCT, by = 'genes'), IR_OV, by = 'genes'),
                                IR_GBM, by = 'genes'),IR_DLBC, by = 'genes'),IR_LUAD,by = 'genes'), IR_BASAL,by = 'genes')
IR_tot<-IR_tot[!duplicated(IR_tot$genes),]
dim(IR_tot)
head(IR_tot)
IR_IR<-merge(IR_tot, IR_L1_Alu, by.x = 'genes', by.y = 'V5')
IR_IR<-IR_IR[!duplicated(IR_IR$genes),c(1,9:13)]
dim(IR_IR)
head(IR_IR)
write.table(IR_IR,'~/Downloads/MIDDLE_tot.txt', sep = '\t')

###################################
require("reticulate")

source_python("pickle_reader.py")
pickle_data<-read_pickle
#####################################
#bedtools intersect -a LAST_LINE.bed -b LAST_SINE.bed -wa|awk '!seen[$0]++'|sort -k1,1 -k2,2n> LAST_L_S_INE.bed
######################################
#install.packages("VennDiagram")
library("VennDiagram")

grid.newpage()
draw.pairwise.venn(area1=64,
                   area2=72,
                   cross.area=38)
######################






















sum_lib<-c(sum(hnM_tot$HTE1)/1000000, sum(hnM_tot$HTE2)/1000000, sum(hnM_tot$HTE3)/1000000, sum(hnM_tot$HTE4)/1000000,
           sum(hnM_tot$HH_5)/1000000, sum(hnM_tot$HH_6)/1000000, sum(hnM_tot$HH5)/1000000, sum(hnM_tot$HH6)/1000000)
for (i in 3:ncol(hnM_tot)){
  hnM_tot[,i]<-hnM_tot[,i]/sum_lib[i-2]
}

for (k in 1:nrow(hnM_tot)){
  hnM_tot[k,3:ncol(hnM_tot)]<-(hnM_tot[k,3:ncol(hnM_tot)]*1000)/hnM_tot$len[k]
}

### 
# convert Ensembl gene id to entrez ID
#library(biomaRt)
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),
  values=hnM_tot$gene,
  mart=mart)
head(genes)
dim(genes)
symbol_HMLE<-merge(genes, hnM_tot, by.x="ensembl_gene_id", by.y="gene")
symbol_HMLE<-symbol_HMLE[!duplicated(symbol_HMLE$ensembl_gene_id),]
head(symbol_HMLE)
# 55259; 20 
dim(symbol_HMLE)
#  cat hallmark_interferon_alpha_response_gene_set.txt hallmark_interferon_gamma_response_gene_set.txt A.J.Scadden_nature_struct_bio.interferon_I_list|awk '!seen[$0]++' > ISG.txt
ISG<-read.table('~/Downloads/ISG.txt')
head(ISG)
hnM_ISG<-merge(symbol_HMLE, ISG, by.x='hgnc_symbol', by.y='V1')
hnM_ISG<-hnM_ISG[!duplicated(hnM_ISG$hgnc_symbol),]
head(hnM_ISG)
dim(hnM_ISG)
ISG_input<-hnM_ISG[,-c(1:4)]
rownames(ISG_input)<-hnM_ISG$hgnc_symbol
################################################################################
################################################################################
################################################################################
FC_input<-ISG_input[,1:4]
colnames(FC_input)<-c('pair1','pair2','pair3', 'pair4')
head(FC_input)
# 5829; 9
dim(FC_input)
for (k in 1:nrow(FC_input)){
  FC_input[k,1]<-(ISG_input[k,2]+0.01)/(ISG_input[k,1]+0.01)
  FC_input[k,2]<-(ISG_input[k,4]+0.01)/(ISG_input[k,3]+0.01)
  FC_input[k,3]<-(ISG_input[k,6]+0.01)/(ISG_input[k,5]+0.01)
  FC_input[k,4]<-(ISG_input[k,8]+0.01)/(ISG_input[k,7]+0.01)
}
mean_KD<-c()
for (i in 1:nrow(FC_input)){
  mean_KD<-c(mean_KD, mean(as.numeric(FC_input[i,1:4])))
}
FC_input$mean_KD<-mean_KD
head(FC_input)
summary(FC_input)
new_rev_trim<-FC_input[rev(order(FC_input$mean_KD)),-ncol(FC_input)]
new_rev_trim<-na.omit(new_rev_trim)
new_rev_trim[new_rev_trim>=2.5]<-2.5
summary(new_rev_trim)
head(new_rev_trim)
rev_rev_trim<-new_rev_trim[rev(order(rowSums(new_rev_trim))),]
### ??? What is the method they use in Fig6C ??? 
library(pheatmap)
library(RColorBrewer)
#install.packages('viridis')
library(viridis)
#breaksList = seq(0, 100, by = 1)
pheatmap(mat=rev_rev_trim[1:233,],
         scale='none',
         #clustering_distance_rows = "euclidean",
         #clustering_method = "complete",
         #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         color=colorRampPalette(c("royalblue","white","red","red3"))(n=30),
         cluster_rows = FALSE,
         cluster_cols = TRUE)















#####################################
hnM_tot<-merge(merge(merge(merge(merge(merge(merge(merge(convert[,2:3], HepG2[,1:2], by='gene'),HNRNPM[,1:2], by='gene'), MATR3[,1:2], by='gene'), SUGP2[,1:2], by='gene'), RBFOX2[,1:2], by='gene'),
                           HNRNPL[,1:2], by='gene'), GTF2F1[,1:2], by='gene'), AKAP8[,1:2], by='gene')
hnM_tot<-hnM_tot[!duplicated(hnM_tot$gene),]
dim(hnM_tot)
head(hnM_tot)
sum_lib<-c(sum(hnM_tot$HepG2)/1000000, sum(hnM_tot$HNRNPM)/1000000, sum(hnM_tot$MATR3)/1000000, sum(hnM_tot$SUGP2)/1000000,
           sum(hnM_tot$RBFOX2)/1000000, sum(hnM_tot$HNRNPL)/1000000, sum(hnM_tot$GTF2F1)/1000000, sum(hnM_tot$AKAP8)/1000000)
for (i in 3:ncol(hnM_tot)){
  hnM_tot[,i]<-hnM_tot[,i]/sum_lib[i-2]
}
#hnM_tot[,9]<-hnM_tot[,9]/sum_lib[8]
hnM_tot_mod<-hnM_tot
for (k in 1:nrow(hnM_tot)){
  hnM_tot_mod[k,3:ncol(hnM_tot_mod)]<-as.numeric((hnM_tot[k,3:ncol(hnM_tot)]*1000))/hnM_tot$len[k]
}
hnM_tot_del<-hnM_tot_mod[rowSums(hnM_tot_mod[,-c(1:2)]>5)>=1,]
# 5998
dim(hnM_tot_del)
head(hnM_tot_del)

####################################
### 
# convert Ensembl gene id to entrez ID
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),
  values=hnM_tot_del$gene,
  mart=mart)
head(genes)
dim(genes)
symbol_HMLE<-merge(genes, hnM_tot_del, by.x="ensembl_gene_id", by.y="gene")
symbol_HMLE<-symbol_HMLE[!duplicated(symbol_HMLE$ensembl_gene_id),]
head(symbol_HMLE)
write.table(symbol_HMLE$hgnc_symbol,'~/Downloads/input_background_file_RBPs.txt', row.names = F, col.names = F, quote = F)
write.table(symbol_HMLE[,c(3,5)],'~/Downloads/rpkm_file.txt', row.names = F, col.names = F, quote = F)

# 55259; 20 
dim(symbol_HMLE)
hnM_ISG<-merge(symbol_HMLE, ISG, by.x='hgnc_symbol', by.y='sample')
hnM_ISG<-hnM_ISG[!duplicated(hnM_ISG$hgnc_symbol),]
head(hnM_ISG)
dim(hnM_ISG)
ISG_input<-hnM_ISG[,-c(1,4,ncol(hnM_ISG))]
symbol_HMLE<-symbol_HMLE[!duplicated(symbol_HMLE$hgnc_symbol),]
ISG_input<-symbol_HMLE[,-c(3:4)]
rownames(ISG_input)<-symbol_HMLE$hgnc_symbol
rownames(ISG_input)<-hnM_ISG$hgnc_symbol
head(ISG_input)
tmp_ISG<-data.frame(cbind('gene'=rownames(ISG_input), ISG_input))
head(tmp_ISG)
write.table(ISG_input,'~/Downloads/ISG_input.txt', sep='\t', quote = F, row.names = F)
