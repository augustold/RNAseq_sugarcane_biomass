
#########################################################################################
## DESeq
#########################################################################################

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("ashr")

library("DESeq2")
library("apeglm")
library("ashr")
library("tximport")
library("stringr")
library("openxlsx")

#Run DESeq2

tx2gene = read.csv("tx2gene.csv", header = T)

contrast_files = list.files("./coldata_contrasts")
contrast = as.data.frame(str_split_fixed(contrast_files, "coldata_", 2))
contrast = as.data.frame(str_split_fixed(contrast$V2, ".txt", 2))
contrast = as.character(contrast$V1)

DESeq_results = sapply(contrast,function(x) NULL)

for (i in contrast) {
  coldata = read.table(paste("./coldata_contrasts/coldata_",i,".txt", sep=""), header = T)
  samples = as.character(coldata$Sample)
  files = file.path("/Users/augustodiniz/Documents/Documentos_CSHL/scga7_and_spon_quants", paste(samples, "_quant", sep=""), "quant.sf")
  names(files) = samples
  all(file.exists(files))
  
  txi = tximport(files, type = "salmon", tx2gene = tx2gene)
  
  sampleTable = data.frame(Class = coldata$Class)
  rownames(sampleTable) = colnames(txi$counts)
  
  dds = DESeqDataSetFromTximport(txi, sampleTable, ~Class)
  
  #pre-filtering
  keep = rowSums(counts(dds)) >= 10
  dds = dds[keep,]
  
  #DEGs
  
  dds = DESeq(dds)
  resultsNames(dds)
  resLFC = lfcShrink(dds, coef=2, type="apeglm")
  resLFC_df = as.data.frame(resLFC)
  resLFC_df$abs_val = abs(resLFC_df$log2FoldChange)
  resLFC_df$genemodel = rownames(resLFC_df)
  resLFC_df = subset(resLFC_df, abs_val >= 2 & padj <= 0.05)
  
  DESeq_results[[i]] = resLFC_df
}

save(DESeq_results, file = "DESeq_results.RData")
write.xlsx(DESeq_results, file = "DEGs_tissue_within_collection_pairs.xlsx", row.names = F, dec = ",")

#=====================================================================================
#
#  UP and DOWN summary
#
#=====================================================================================

library(dplyr)
library(stringr)
require(openxlsx)
load("DESeq_results.RData")

#UP scga7
scga7_UP_DEGs = vector("list", length(DESeq_results))
names(scga7_UP_DEGs) = paste("UP", names(DESeq_results))
for (i in 1:(length(scga7_UP_DEGs))) {
  x = as.data.frame(DESeq_results[i])
  y = as.data.frame(str_split_fixed(x[,7], "_", 2))
  row.names(y) = row.names(x)
  y = subset(y, V1 == "scga7")
  y = as.character(row.names(y))
  x = x[y,]
  names(x) = c("baseMean", "log2FoldChange", "lfcSE", "pvalue",
               "padj", "abs_val", "genemodel")
  UP = subset(x, log2FoldChange > 0)
  scga7_UP_DEGs[[i]] = UP
}

#DOWN scga7
scga7_DOWN_DEGs = vector("list", length(DESeq_results))
names(scga7_DOWN_DEGs) = paste("DOWN", names(DESeq_results))
for (i in 1:(length(scga7_DOWN_DEGs))) {
  x = as.data.frame(DESeq_results[i])
  y = as.data.frame(str_split_fixed(x[,7], "_", 2))
  row.names(y) = row.names(x)
  y = subset(y, V1 == "scga7")
  y = as.character(row.names(y))
  x = x[y,]
  names(x) = c("baseMean", "log2FoldChange", "lfcSE", "pvalue",
               "padj", "abs_val", "genemodel")
  DOWN = subset(x, log2FoldChange < 0)
  scga7_DOWN_DEGs[[i]] = DOWN
}

#UP spon
spon_UP_DEGs = vector("list", length(DESeq_results))
names(spon_UP_DEGs) = paste("UP", names(DESeq_results))
for (i in 1:(length(spon_UP_DEGs))) {
  x = as.data.frame(DESeq_results[i])
  y = as.data.frame(str_split_fixed(x[,7], "_", 2))
  row.names(y) = row.names(x)
  y = subset(y, V1 != "scga7")
  y = as.character(row.names(y))
  x = x[y,]
  names(x) = c("baseMean", "log2FoldChange", "lfcSE", "pvalue",
               "padj", "abs_val", "genemodel")
  UP = subset(x, log2FoldChange > 0)
  spon_UP_DEGs[[i]] = UP
}

#DOWN spon
spon_DOWN_DEGs = vector("list", length(DESeq_results))
names(spon_DOWN_DEGs) = paste("DOWN", names(DESeq_results))
for (i in 1:(length(spon_DOWN_DEGs))) {
  x = as.data.frame(DESeq_results[i])
  y = as.data.frame(str_split_fixed(x[,7], "_", 2))
  row.names(y) = row.names(x)
  y = subset(y, V1 != "scga7")
  y = as.character(row.names(y))
  x = x[y,]
  names(x) = c("baseMean", "log2FoldChange", "lfcSE", "pvalue",
               "padj", "abs_val", "genemodel")
  DOWN = subset(x, log2FoldChange < 0)
  spon_DOWN_DEGs[[i]] = DOWN
}

scga7_N_UP_DEGs = vector()
scga7_N_DOWN_DEGs = vector()

spon_N_UP_DEGs = vector()
spon_N_DOWN_DEGs = vector()

for (i in 1:8) {
  scga7_N_UP_DEGs[i] = as.numeric(nrow(as.data.frame(scga7_UP_DEGs[[i]])))
  scga7_N_DOWN_DEGs[i] = as.numeric(nrow(as.data.frame(scga7_DOWN_DEGs[[i]])))
  spon_N_UP_DEGs[i] = as.numeric(nrow(as.data.frame(spon_UP_DEGs[[i]])))
  spon_N_DOWN_DEGs[i] = as.numeric(nrow(as.data.frame(spon_DOWN_DEGs[[i]])))
}

DEG_summary = as.data.frame(cbind(scga7_N_UP_DEGs,
                                  scga7_N_DOWN_DEGs,
                                  spon_N_UP_DEGs,
                                  spon_N_DOWN_DEGs))
row.names(DEG_summary) = names(DESeq_results)


#=====================================================================================
#
#  DEG annotation
#
#=====================================================================================

library(dplyr)
library(stringr)
require(openxlsx)
load("DESeq_results.RData")

ref_track = read.csv("genome_ref_track.csv", header = T)
row.names(ref_track) = ref_track[,1]
scga7 = subset(ref_track, REF == "scga7")
spon = subset(ref_track, REF == "spon")

annot_scga7 = read.delim("KEGG_annot_scga7_v2.txt", header = F)
names(annot_scga7) = c("GENEID", "KOID", "GENESYMBOL", "DESCRIPTION", "EC")
row.names(annot_scga7) = annot_scga7[,1]
annot_scga7 = inner_join(scga7, annot_scga7, "GENEID")
TF_scga7 = read.delim("TF_scga7.txt", header = F)
names(TF_scga7) = c("GENEID", "TF")
annot_scga7 = full_join(annot_scga7, TF_scga7,  "GENEID")

annot_spon = read.delim("KEGG_annot_spont_v2.txt", header = F)
names(annot_spon) = c("GENEID", "KOID", "GENESYMBOL", "DESCRIPTION", "EC")
row.names(annot_spon) = annot_spon[,1]
annot_spon = inner_join(spon, annot_spon, "GENEID")
TF_spon = read.delim("TF_spon.txt", header = F)
names(TF_spon) = c("GENEID", "TF")
annot_spon = full_join(annot_spon, TF_spon,  "GENEID")

#Annotated DEGs scga7
scga7_annot_DEGs = vector("list", length(DESeq_results))
names(scga7_annot_DEGs) = names(DESeq_results)
for (i in 1:(length(scga7_annot_DEGs))) {
  x = as.data.frame(DESeq_results[i])
  y = as.data.frame(str_split_fixed(x[,7], "_", 2))
  row.names(y) = row.names(x)
  y = subset(y, V1 == "scga7" | V1 == "mikado" | V1 == "novel")
  y = as.character(row.names(y))
  x = x[y,]
  names(x) = c("baseMean", "log2FoldChange", "lfcSE", "pvalue",
               "padj", "abs_val", "GENEID")
  w = full_join(x, annot_scga7, "GENEID")
  w = w[,-c(1:6)]
  x = inner_join(x, w, "GENEID")
  scga7_annot_DEGs[[i]] = x
}

#Annotated DEGs spon
spon_annot_DEGs = vector("list", length(DESeq_results))
names(spon_annot_DEGs) = names(DESeq_results)
for (i in 1:(length(spon_annot_DEGs))) {
  x = as.data.frame(DESeq_results[i])
  y = as.data.frame(str_split_fixed(x[,7], "_", 2))
  row.names(y) = row.names(x)
  y = subset(y,  V1 != "scga7")
  y = subset(y,  V1 != "mikado")
  y = subset(y,  V1 != "novel")
  y = as.character(row.names(y))
  x = x[y,]
  names(x) = c("baseMean", "log2FoldChange", "lfcSE", "pvalue",
               "padj", "abs_val", "GENEID")
  w = full_join(x, annot_spon, "GENEID")
  w = w[,-c(1:6)]
  x = inner_join(x, w, "GENEID")
  spon_annot_DEGs[[i]] = x
}

save(scga7_annot_DEGs, file = "scga7_annot_DEGs.RData")
save(spon_annot_DEGs, file = "spon_annot_DEGs.RData")

write.xlsx(scga7_annot_DEGs, file = "DEGs_scga7_annot.xlsx", row.names = F, dec = ",")
write.xlsx(spon_annot_DEGs, file = "DEGs_spon_annot.xlsx", row.names = F, dec = ",")


#=====================================================================================
#
#  GSEA
#
#=====================================================================================

library(fgsea)
require(openxlsx)
library(data.table)
library(ggplot2)
library(BiocParallel)
register(SerialParam())

#scga7

load("scga7_annot_DEGs.RData")
contrasts = names(scga7_annot_DEGs)

gmt_scga7 = read.delim("gmt_gramene_scga7.txt", header = T)
gmt_scga7 = unique(gmt_scga7)
gmt_scga7 = na.omit(gmt_scga7)
Pathways = as.character(unique(gmt_scga7$term))

Pathways_list_scga7 = sapply(Pathways,function(x) NULL)

for (i in Pathways) {
  x = subset(gmt_scga7, term == i)
  x = as.character(x$gene)
  Pathways_list_scga7[[i]] = x
}

gsea_results_scga7 = sapply(contrasts,function(x) NULL)

for (i in contrasts) {
  x = as.data.frame(scga7_annot_DEGs[[i]])
  ranks = x$log2FoldChange
  names(ranks) = x$GENEID
  
  fgseaRes = fgsea(pathways = Pathways_list_scga7, 
                   stats    = ranks,
                   maxSize  = 2000)
  fgseaRes = subset(fgseaRes, pval <= 0.05)
  gsea_results_scga7[[i]] = fgseaRes[order(pval), ]
}

save(gsea_results_scga7, file = "gsea_results_scga7.RData")

#spon

load("spon_annot_DEGs.RData")
contrasts = names(spon_annot_DEGs)

gmt_spon = read.delim("gmt_gramene_spon.txt", header = T)
gmt_spon = unique(gmt_spon)
gmt_spon = na.omit(gmt_spon)
Pathways = as.character(unique(gmt_spon$term))

Pathways_list_spon = sapply(Pathways,function(x) NULL)

for (i in Pathways) {
  x = subset(gmt_spon, term == i)
  x = as.character(x$gene)
  Pathways_list_spon[[i]] = x
}

gsea_results_spon = sapply(contrasts,function(x) NULL)

for (i in contrasts) {
  x = as.data.frame(spon_annot_DEGs[[i]])
  ranks = x$log2FoldChange
  names(ranks) = x$GENEID
  
  fgseaRes = fgsea(pathways = Pathways_list_spon, 
                   stats    = ranks,
                   maxSize  = 2000)
  fgseaRes = subset(fgseaRes, pval <= 0.05)
  gsea_results_spon[[i]] = fgseaRes[order(pval), ]
}

save(gsea_results_spon, file = "gsea_results_spon.RData")

write.xlsx(gsea_results_scga7, file = "gsea_results_scga7.xlsx", row.names = F, dec = ",")
write.xlsx(gsea_results_spon, file = "gsea_results_spon.xlsx", row.names = F, dec = ",")


#=====================================================================================
#
#  DEGs Database
#
#=====================================================================================

library(dplyr)
library(stringr)

#DEG List

load("DESeq_results.RData")

DEGs_database = data.frame()

contrasts = as.character(names(DESeq_results))

for (i in contrasts) {
  
  x = as.data.frame(DESeq_results[i])
  x = as.data.frame(x[,7])
  names(x) = "GENEID"
  
  temp_name = as.data.frame(str_split_fixed(i, "_", 2))
  
  if (temp_name[1,1] == "Tissue") {
    
    temp_name2 = as.data.frame(str_split_fixed(i, "_", 2))
    x$ContrastGroup = rep(as.character(temp_name2[1,1]), nrow(x))
    x$Contrast = rep(as.character(temp_name2[1,2]), nrow(x))
    x$Tissue = rep(NA, nrow(x))
    x$Collection = rep(NA, nrow(x))
    
  }
  
  else {
    
    temp_name3 = as.data.frame(str_split_fixed(i, "_", 4))
    x$ContrastGroup = rep(as.character(temp_name3[1,1]), nrow(x))
    x$Contrast = rep(as.character(temp_name3[1,4]), nrow(x))
    x$Tissue = rep(as.character(temp_name3[1,2]), nrow(x))
    x$Collection = rep(as.character(temp_name3[1,3]), nrow(x))
    
  }
  
  DEGs_database = rbind(DEGs_database, x)
  
}

save(DEGs_database, file = "DEGs_database.RData")

# Annotation

ref_track = read.csv("genome_ref_track.csv", header = T)
row.names(ref_track) = ref_track[,1]
scga7 = subset(ref_track, REF == "scga7")
spon = subset(ref_track, REF == "spon")

annot_scga7 = read.delim("KEGG_annot_scga7_v2.txt", header = F)
names(annot_scga7) = c("GENEID", "KOID", "GENESYMBOL", "DESCRIPTION", "EC")
row.names(annot_scga7) = annot_scga7[,1]
annot_scga7 = inner_join(scga7, annot_scga7, "GENEID")
TF_scga7 = read.delim("TF_scga7.txt", header = F)
names(TF_scga7) = c("GENEID", "TF")
annot_scga7 = full_join(annot_scga7, TF_scga7,  "GENEID")

annot_spon = read.delim("KEGG_annot_spont_v2.txt", header = F)
names(annot_spon) = c("GENEID", "KOID", "GENESYMBOL", "DESCRIPTION", "EC")
row.names(annot_spon) = annot_spon[,1]
annot_spon = inner_join(spon, annot_spon, "GENEID")
TF_spon = read.delim("TF_spon.txt", header = F)
names(TF_spon) = c("GENEID", "TF")
annot_spon = full_join(annot_spon, TF_spon,  "GENEID")

Functional_annotation = rbind(annot_scga7, annot_spon)
save(Functional_annotation, file = "Functional_annotation.RData")

# Subset genes of interest

MIG = subset(annot_scga7, DESCRIPTION == "sucrose synthase")
MIG = as.character(MIG$GENEID)

MIG_DEG = DEGs_database[(DEGs_database$GENEID) %in% MIG,]



#=====================================================================================
#
#  Heatmap
#
#=====================================================================================
library("pheatmap")
library("DESeq2")


load("DEGs_database.RData")
load("Functional_annotation.RData")


load("contrasts_1.RData")
load("vsd.RData")
exp_vsd = as.data.frame(assay(vsd))
rm(vsd)
exp_vsd = exp_vsd[,as.character(sample_annot$Sample)]

varlist = as.character(names(contrasts_1))

L1_C1_Low_vs_High = read.table("coldata_L1_C1.txt", header = T)
L1_C2_Low_vs_High = read.table("coldata_L1_C2.txt", header = T)
L1_C3_Low_vs_High = read.table("coldata_L1_C3.txt", header = T)
I1_C1_Low_vs_High = read.table("coldata_I1_C1.txt", header = T)
I1_C2_Low_vs_High = read.table("coldata_I1_C2.txt", header = T)
I1_C3_Low_vs_High = read.table("coldata_I1_C3.txt", header = T)
I5_C2_Low_vs_High = read.table("coldata_I5_C2.txt", header = T)
I5_C3_Low_vs_High = read.table("coldata_I5_C3.txt", header = T)

for (i in varlist) {
  select = as.data.frame(contrasts_1[i])
  select = as.character(select[,ncol(select)])
  
  df = get(i)
  row.names(df) = df$Sample
  df = df[,-c(1,4,6,8)]
  
  mat_exp = exp_vsd[select,as.character(row.names(df))]
  thresh = mat_exp  > 0
  keep = rowSums(thresh) >= 1
  exp_vsd_keep = mat_exp[keep,]
  exp_vsd_keep = na.omit(exp_vsd_keep)
  exp_vsd_keep[exp_vsd_keep == 0] = 0.001
  
  png(file = paste(i, ".png", sep=""))
  pheatmap(as.matrix(log10(exp_vsd_keep)),
           cluster_rows=T, cluster_cols=T,
           show_rownames=FALSE,
           show_colnames=FALSE,
           annotation_col=df)
  dev.off()
}









