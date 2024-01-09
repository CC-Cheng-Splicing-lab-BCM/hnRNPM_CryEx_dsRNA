######low hnRNPM associates w. immune responses upregulation among cancers######
### 1. MsigDB hallmark (hnM low vs. high)
library(fgsea)
library(ggplot2)
# hallmark gmt txt file downloadable through GSEA website
filename = "~/Downloads/h.all.v7.4.symbols.gmt.txt"
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
###########################Using TCGA BRCA as an example######################################################
cancer_types = c('BRCA')
# define working directory
base_dir = '~/Downloads/'

# counts table
counts_filename = sapply(cancer_types, function(x) paste(base_dir,x,'_out_counts.txt',sep = ''))                                 
# counts metadata
counts_meta = sapply(cancer_types, function(x) paste('~/Downloads/gdc_', x, '_counts.tsv', sep = ''))
###########################
i=1
BRCA_meta_ori<-read.table(counts_meta[[i]], sep = '\t', header = T)
colnames(BRCA_meta_ori)<-c('File.ID', 'File.Name', 'Data.Category', 'Data.Type', 'Project.ID', 
                          'Case.ID', 'Sample.ID', 'Sample.Type')
BRCA_meta<-BRCA_meta_ori[BRCA_meta_ori$Sample.Type=="Primary Tumor"|BRCA_meta_ori$Sample.Type=="Primary Blood Derived Cancer - Peripheral Blood",]
BRCA_meta$File.Name<-gsub("-", ".", substring(BRCA_meta$File.Name,1,36))
rownames(BRCA_meta)<-BRCA_meta$File.Name
BRCA_meta$Sample.ID<-substring(BRCA_meta$Sample.ID,1,15)
head(BRCA_meta)
dim(BRCA_meta)
# # ############################
# # to extract BRCA different subtypes samples:
# metadata<-read.csv('~/Downloads/clinical_data_BRCA.txt', sep = '\t')
# # PAM50Call_RNAseq
# unique(metadata$PAM50Call_RNAseq)
# LumA_ori<-metadata[metadata$PAM50Call_RNAseq=="LumB",c('sampleID','sample_type')]
# LumA<-LumA_ori[LumA_ori$sample_type=="Primary Tumor",]
# # get rid of a metastatic sample
# dim(LumA)
# 
# BRCA_meta_new<-merge(BRCA_meta, LumA, by.x="Sample.ID", by.y="sampleID")
# BRCA_meta_new<-BRCA_meta_new[!duplicated(BRCA_meta_new$Sample.ID),]
# dim(BRCA_meta_new)
# head(BRCA_meta_new)
################################################################################
######lcpmtmm normalization
BRCA_count<-read.table(counts_filename[[i]], sep = '\t', header = T)
BRCA_count$Ensembl_gene_id<-substring(BRCA_count$Ensembl_gene_id,1,15)
BRCA_count<-BRCA_count[!duplicated(BRCA_count$Ensembl_gene_id),]
rownames(BRCA_count)<-BRCA_count$Ensembl_gene_id
colnames(BRCA_count)<-gsub("X","", colnames(BRCA_count))
colnames(BRCA_count)<-c('Ensembl_gene_id', gsub(".htseq.counts","",colnames(BRCA_count[,-1])))
colnames(BRCA_count)<-c('Ensembl_gene_id', gsub("\\.","-",colnames(BRCA_count[,-1])))
dim(BRCA_count)
head(BRCA_count[1:3,1:3])
BRCA_meta$File.Name<-gsub('\\.', '-', BRCA_meta$File.Name)
BRCA_count_trim<-BRCA_count[,c(BRCA_meta$File.Name)]
# BRCA_meta_new$File.Name<-gsub('\\.', '-', BRCA_meta_new$File.Name)
# BRCA_count_trim<-BRCA_count[,c(BRCA_meta_new$File.Name)]
head(BRCA_count_trim[,1:5])
dim(BRCA_count_trim)
###
#BiocManager::install('DESeq2')
library('DESeq2')
#BiocManager::install('edgeR')
library('edgeR')

input<-na.omit(BRCA_count_trim)

conditions<-na.omit(factor(colnames(BRCA_count_trim)))

colData<-na.omit(data.frame(sampleNames=colnames(input), conditions=conditions))
row.names(colData)<-colData$sampleNames
dds<-DESeqDataSetFromMatrix(countData = as.matrix(input),
                            colData = colData,
                            design= ~conditions)
y<-DGEList(counts(dds))
y$samples
#### get rid of low abundance genes
keep<-rowSums(cpm(counts(dds))>5/(min(y$samples$lib.size)/10^6))>(length(y$samples$lib.size)/2)
dds_kpt<-dds[keep,]
dim(dds_kpt)
##################
ym<-DGEList(counts(dds_kpt))
ym<-calcNormFactors(ym, method="TMM")
lcpmtmm.rna<-cpm(ym, log=T, normalized.lib.sizes=T)
colnames(lcpmtmm.rna)<-gsub("X","",colnames(lcpmtmm.rna))
head(lcpmtmm.rna[,1:5])
# hnRNPM: ENSG00000099783
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
##################
group <- c(rep("LQ",length(lq)), rep("HQ",length(hq)))
design <- model.matrix(~1+group)
counts<-ym$counts[,c(lq,hq)]
dge<-DGEList(counts = counts)
keep<-filterByExpr(dge, design)
dge_trim<-dge[keep,]
dge_new<-calcNormFactors(dge_trim)
v<-voom(dge_new, design, plot = F)
v.fit<-lmFit(v, design)
v.fit<-eBayes(v.fit)
v.res<-topTable(v.fit, coef = ncol(design), number = nrow(v$E))
v.res_id<-data.frame(cbind("ensembl_gene_id"=substr(rownames(v.res),1,15),v.res))
head(v.res_id)
###convert to gene symbol
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),
  values=v.res_id$ensembl_gene_id,
  mart=mart)
head(genes)
dim(genes)
symbol<-merge(v.res_id, genes, by='ensembl_gene_id')
symbol<-symbol[!duplicated(symbol$hgnc_symbol),]
head(symbol)
dim(symbol)
################################################################################
input<-as.numeric(symbol$logFC)
names(input)<-symbol$hgnc_symbol
head(input)

fgseaRes <- fgsea(pathways=all_pathway, 
                  stats=input,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)
fgseaRes[fgseaRes$padj<0.05,c(1,3,5)]
# Repeat the above for each TCGA cancer and save the hallmark scores in its individual file
write.table(fgseaRes[,c(1,5,2,3)],'~/Downloads/Hallmark_BRCA_Jul_2022.txt', sep='\t', row.names = FALSE)
################################################################################
### plot the bubble plot
library(dplyr)
library(reshape2)
library(data.table)

#########################################################################################
# load files containing go enrichment results (merge all TCGA cancers in one bubble plot)

# header of tables under go_enrich_filenames:
# unique_query_GO_terms	occurrences_in_query	enrichment_p_values	FDR
#install.packages("data.table")
library("data.table")
cancer_types<-c('BRCA_LumA','BRCA_LumB','BRCA_Her2','BRCA_Basal','ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC',
                'ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD',
                'PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')

go_enrich_df_list = list()
base_dir = '~/Downloads/Hallmark_'
# hallmark table
hallmark_filename = sapply(cancer_types, function(x) paste(base_dir, x,
                                                           '_Jul_2022.txt', sep = ''))
for (i in seq(length(cancer_types))){
  enrich_df = fread(hallmark_filename[[i]], sep='\t', header=T, stringsAsFactors = F, check.names = F)
  enrich_df$cancer_type = cancer_types[i]
  go_enrich_df_list[[i]] = enrich_df
  
}

sig_enrich_list = lapply(go_enrich_df_list, function(x) {x = x %>% filter(NES >= -5 & padj <= 1)
x$neglogq = -log10(x$padj)
return(x)})


anno_sig_enrich_list = sig_enrich_list

bp_sig_list = sig_enrich_list
combined_all_bp_enrich_df = do.call(rbind, bp_sig_list)

count_mult_bp_enrich_df = as.data.frame(table(combined_all_bp_enrich_df$pathway)) 
bp_mult_enrich_df = combined_all_bp_enrich_df
# replace BRCAording to enrichment
bp_immune_viral_terms = c("ALLOGRAFT_REJECTION","COAGULATION","COMPLEMENT","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE",
                          "IL6_JAK_STAT3_SIGNALING","INFLAMMATORY_RESPONSE")
bp_cancer_terms = c("APOPTOSIS","HYPOXIA","PROTEIN_SECRETION","UNFOLDED_PROTEIN_RESPONSE","REACTIVE_OXYGEN_SPECIES_PATHWAY",
                    "E2F_TARGETS","G2M_CHECKPOINT","MYC_TARGETS_V1","MYC_TARGETS_V2","P53_PATHWAY","MITOTIC_SPINDLE")
bp_emt_terms = c("ADIPOGENESIS","ANGIOGENESIS","EPITHELIAL_MESENCHYMAL_TRANSITION","MYOGENESIS","SPERMATOGENESIS","PANCREAS_BETA_CELLS")
bp_other_terms = c("APICAL_JUNCTION","APICAL_SURFACE","PEROXISOME","DNA_REPAIR","UV_RESPONSE_DN","UV_RESPONSE_UP","BILE_ACID_METABOLISM",
                   "CHOLESTEROL_HOMEOSTASIS","FATTY_ACID_METABOLISM","GLYCOLYSIS","HEME_METABOLISM","OXIDATIVE_PHOSPHORYLATION","XENOBIOTIC_METABOLISM",
                   "ANDROGEN_RESPONSE","ESTROGEN_RESPONSE_EARLY","ESTROGEN_RESPONSE_LATE","IL2_STAT5_SIGNALING","KRAS_SIGNALING_UP",
                   "KRAS_SIGNALING_DN","MTORC1_SIGNALING","NOTCH_SIGNALING","PI3K_AKT_MTOR_SIGNALING","HEDGEHOG_SIGNALING",
                   "TGF_BETA_SIGNALING","TNFA_SIGNALING_VIA_NFKB","WNT_BETA_CATENIN_SIGNALING")

# groups for coloring text (replace or remove as desired)
bp_terms_grouped_df = data.frame(term = c(bp_immune_viral_terms, bp_cancer_terms, bp_emt_terms, bp_other_terms),
                                 gogroup = c(rep('immune/viral', length(bp_immune_viral_terms)),
                                             rep('cancer', length(bp_cancer_terms)),
                                             rep('EMT', length(bp_emt_terms)),
                                             rep('other', length(bp_other_terms))))
library(RColorBrewer)
text_colors = brewer.pal(length(unique(bp_terms_grouped_df$gogroup))-1, 'Dark2')

bp_terms_grouped_df$text_col = sapply(bp_terms_grouped_df$gogroup, function(x) ifelse(x == "immune/viral",text_colors[1],
                                                                                      ifelse(x == "cancer", text_colors[2],
                                                                                             ifelse(x == "EMT",text_colors[3], 'black'))))

plot_bp_mult_enrich_df = merge(bp_terms_grouped_df, bp_mult_enrich_df, by.x='term', by.y='pathway')
plot_bp_mult_enrich_df$term = factor(plot_bp_mult_enrich_df$term, levels=rev(bp_terms_grouped_df$term))

plot_bp_mult_enrich_df$cancer_type = factor(plot_bp_mult_enrich_df$cancer_type, levels = cancer_types,
                                            labels = sapply(cancer_types, toupper))
plot_bp_mult_enrich_df = plot_bp_mult_enrich_df[order(plot_bp_mult_enrich_df$term),]
tmp_df<-plot_bp_mult_enrich_df$neglogq
tmp_df[tmp_df< -log10(0.05)]<-0
plot_bp_mult_enrich_df_new<-data.frame(plot_bp_mult_enrich_df[,-ncol(plot_bp_mult_enrich_df)], 'neglogq'=tmp_df)
head(plot_bp_mult_enrich_df_new)
dim(plot_bp_mult_enrich_df_new)
dim(plot_bp_mult_enrich_df)
# make bubble plot
library(ggplot2)
p_allDE = ggplot(plot_bp_mult_enrich_df_new, aes(x=term, y=cancer_type, size=neglogq,colour=NES)) + 
  geom_point() +#, shape=21) + #
  scale_color_gradient2(low='#0000FF', mid='#ccccccff',high='#FF0000')+
  xlab('') + ylab('') + coord_flip() + scale_size(name='-log10(q-value)') + ##scale_fill_brewer(palette='Dark2', name='') + #
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), #,
        plot.title=element_text(size=12, face="bold"),
        axis.text.x = element_text(angle = 30, vjust=1, hjust=1),
        axis.text.y = element_text(color=rev(bp_terms_grouped_df$text_col)), #order_plot_mult_edit_enrich_df$text_col),
        legend.text=element_text(size=12), legend.key=element_rect(fill="white"), #legend.position = "none",
        legend.title=element_text(size=12), #,face="bold"),
        strip.text.x = element_text(size=12), #axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray"), #panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"), #
        axis.line = element_line(colour = "black"))


ggsave(plot=p_allDE, width=12, height=11, dpi=300, 
       filename = paste('~/Downloads/', 'dsRNA_CryEx_hallmark_Jul_2022.pdf', sep=''),
       useDingbats=F)


