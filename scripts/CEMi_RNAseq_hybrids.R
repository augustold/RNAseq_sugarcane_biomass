
#=====================================================================================
#
#  CEMiTool
#
#=====================================================================================

#Necessary only when running on helix
#.libPaths( c( "/home/projects/R/x86_64-pc-linux-gnu-library/4.0" , .libPaths() ) )

setwd("/Users/augustodiniz/Documents/Documentos_CSHL/RNAseq_Analysis/Hybrids/Salmon_Final_Mikado_V2_2020/CEMi_biomass")

library("CEMiTool")
library("ggplot2")
library("DESeq2")
library("WGCNA")

tissue = c("L1", "I1", "I5")
genome = c("scga7", "spon")

#Load expression data
load("dds.RData")
exp = as.data.frame(assay(dds))

#exp[exp <= 0] = NA

#Load genome reference track file
genome_ref_track = read.csv("genome_ref_track.csv", header = T)
row.names(genome_ref_track) = genome_ref_track[,1]

#Load sample annotation
sample_annot = read.csv("sample_annot.csv", header = T)
sample_annot = subset(sample_annot, Type == "Hybrid")

##For testing
#i = "L1"
#j = "spon"

for (i in tissue) {
  for (j in genome) {
    #Set root directory
    setwd("/Users/augustodiniz/Documents/Documentos_CSHL/RNAseq_Analysis/Hybrids/Salmon_Final_Mikado_V2_2020/CEMi_biomass")
    
    #Set up expression matrix
    select = subset(genome_ref_track, REF == j)
    
    sample_annot_temp = subset(sample_annot, Exp == "Energy cane" & Tissue == i)
    sample_annot_temp = subset(sample_annot_temp,
                               ID == "C411" | 
                                 ID == "C418" | 
                                 ID == "C419" | 
                                 ID == "C429" )
    row.names(sample_annot_temp) = sample_annot_temp$Sample
    
    #Get expression matrix
    exp_temp = exp[as.character(row.names(select)),as.character(sample_annot_temp$Sample)]
    exp_temp[exp_temp <= 0] = NA
    
    #Filter expression matrix
    exp_matrix = t(as.matrix(exp_temp))
    gsg = goodSamplesGenes(exp_matrix, verbose = 3, minFraction = 9/10);
    gsg$allOK
    datExpr_A = exp_matrix[, gsg$goodGenes]
    dim(datExpr_A)
    datExpr_B = as.data.frame(t(datExpr_A))
    datExpr_B[is.na(datExpr_B)] = 0
    
    #Set up annotation info
    SampleName = as.character(sample_annot_temp$Sample)
    Class = as.character(sample_annot_temp$Class3)
    annot = as.data.frame(cbind(SampleName, Class))
    gmt = read.delim(paste("/Users/augustodiniz/Documents/Documentos_CSHL/RNAseq_Analysis/Hybrids/Salmon_Final_Mikado_V2_2020/CEMi_biomass/gmt_gramene_",j,".txt", sep = ""), header = T)
    gmt = unique(gmt)
    
    #Set up directory
    setwd(paste("./",i,"/",j,sep = ""))
    
    #Run cemitool
    cem <- cemitool(datExpr_B, annot, gmt,
                    apply_vst = TRUE,
                    force_beta = TRUE,
                    verbose=TRUE)
    
    dev.set(dev.next())
    
    # create report as pdf and html documents
    generate_report(cem, directory="./Report", force=TRUE,
                    output_format=c( "html_document"))
    
    # write analysis results into files
    write_files(cem, directory="./Tables", force=TRUE)
    
    # save all plots
    save_plots(cem, "all", directory="./Plots", force=TRUE)
    
    diagnostic_report(cem, directory="./Diagnostics", force=TRUE)
    
    save(cem, file = "cem.RData")
    
  }
}

rm(list = ls())
