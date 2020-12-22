##Sample info
#sample_annot = read.csv("sample_annot.csv", header = T)
#names(sample_annot)
#sample_annot = subset(sample_annot, Exp == "Energy cane")
#coldata = sample_annot[,-c(1,3)]

#write.table(coldata, file = "coldata.txt", row.names = F)
coldata = read.table("coldata.txt", header = T)
row.names(coldata) = coldata[,1]
levels(factor(coldata$ID))

tx2gene = read.csv("tx2gene.csv", header = T)

library("DESeq2")
library("apeglm")
library("ashr")
library("tximport")
library("vsn")
library("ggplot2")

samples = as.character(coldata$Sample)
files = file.path("/Users/augustodiniz/Documents/Documentos_CSHL/scga7_and_spon_quants", paste(samples, "_quant", sep=""), "quant.sf")
names(files) = samples
all(file.exists(files))

txi = tximport(files, type = "salmon", tx2gene = tx2gene)

dds = DESeqDataSetFromTximport(txi,
                                coldata,
                                ~ ID + Tissue + Collection)

#pre-filtering
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]
dds

save(dds, file = "dds.RData")

##Variance stabilizing transformation
vsd = vst(dds, blind=FALSE)
exp_vsd = as.data.frame(assay(vsd))
#meanSdPlot(assay(vsd))

save(vsd, file = "vsd.RData")
#write.table(exp_vsd, file = "exp_vsd.txt")

load("vsd.RData")

#PCA
pcaData = plotPCA(vsd, intgroup=c("ID", "Collection"), returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))
pcaData = cbind(pcaData, coldata$Class2, coldata$Tissue)
names(pcaData) = c("PC1","PC2","group","ID","Collection","name","Class2","Tissue")
ggplot(pcaData, aes(PC1, PC2, color=Tissue)) +
  geom_point(size=3, alpha = 0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(breaks = c("I1", "I5", "L1"), 
                      values=c("#00BFC4", "#F8766D", "#7CAE00"))

pdf("PCA_all_samples.pdf", width = 4, height = 4)
ggplot(pcaData, aes(PC1, PC2, color=Tissue)) +
  geom_point(size=3, alpha = 0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(breaks = c("I1", "I5", "L1"), 
                     values=c("#00BFC4", "#F8766D", "#7CAE00"))
dev.off()


library("RColorBrewer")
library("pheatmap")

#Expression Heatmap

select = order(rowMeans(counts(dds,normalized=FALSE)),
                decreasing=TRUE)[1:200]
df = coldata[,c(2,3,5,6)]

pdf("Heatmap_all_samples.pdf", width = 5, height = 7)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F, show_colnames=F,
         cluster_cols=T, annotation_col=df, cutree_rows = 6,cutree_cols = 3)
dev.off()

#Sample dist

sampleDists = dist(t(assay(vsd)))

sampleDistMatrix = as.matrix(sampleDists)
#rownames(sampleDistMatrix) = paste(vsd$ID, vsd$Tissue, vsd$Collection sep="-")
#colnames(sampleDistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
