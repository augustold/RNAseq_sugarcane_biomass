## Overview:

### Pipeline

![Bioinformatics pipeline](https://github.com/augustold/RNAseq_sugarcane_biomass/blob/main/Bioinformatics_pipeline.png)

### This repository contains the scripts for:

**A.** Update genome annotation using the MIKADO pipeline

**1.** Gene prediction
- Generating input data required for running MIKADO pipeline
- Running MIKADO pipeline
- Running PASA pipeline

**2.** Functional annotation
- PlantTFDB
- BlastKOALA
- GRAMENE

**B.** RNAseq anlaysis;

**1.** Exploratory analysis

**2.** Differential expression analysis

**3.** Co-expression analysis


### In order to do that, scripts must be used in the following order:

1. Build genome index for 1st-pass: **star_1stpass.sh**
2. Align reads: **star_align_1stpass.sh**
3. Build genome index for 2nd-pass: **star_2ndpass.sh**
4. Realign reads using new genome index: **star_align_2ndpass.sh**
5. Merge bam files: **merge_bam.sh**
6. Build reference-guided transcriptomes:
- Stringtie: **stringtie.sh**
- Psiclass: **psiclass.sh**
- Cufflinks: **cufflinks.sh**
- Trinity: **trinity.sh**
- Align Trinity transcripts back to the genome: **GMAP_build.sh** and **GMAP_run.sh**
7. Porticullis: *run this step on brie.cshl.edu*
8. Mikado: **mikado.sh**
9. BRAKER: **braker.sh**
10. PASA: **seqclean.sh**, **alignAssembly.sh**, **annotLoadandCompare.sh** and **abinitiotraining.sh**
11. RNAseq exploratory: **Exploratory_script_all_data.R
12. RNAseq DEG: **DEGs.R**
13. RNAseq Co-expression: **CEMi_RNAseq_hybrids.R**

Follow the steps described in 'Workflow.md' file.

### Genome reference

Assembly of the 373k gene space of the polyploid sugarcane genom: DOI https://doi.org/10.1093/gigascience/giz129

Genome data: http://gigadb.org/dataset/100655

### RNA seq data

Samples were collected from a Brazilian sugarcane hybrid variety - SP80-3280 - grown under field conditions at 12 months after planting. Collected tissues include leaf +1 (L1), upper internode (I1), young internode (I5) and mature internode (I9) with three biological replicates per tissue.

Stranded libraries were prepared and sequenced on Illumina Hiseq platform to obtain 150 bp paired-end reads. Each sample holds ~ 20M reads.

### Script run time

**Script**|**threads**|**m\_mem\_free**|**Run Time**|**Days**
:-----:|:-----:|:-----:|:-----:|:-----:
star\_1stpass.sh|16|5|01:08:57| 
star\_align\_1stpass.sh|16|5|22:14:52|1
star\_2ndpass.sh|16|5|01:58:04| 
star\_align\_2ndpass.sh|16|5|21:06:33|1
merge\_bam.sh|5|5|23:07:49|1
stringtie.sh|16|5|72:12:04|3
psiclass.sh|16|5|23:58:41| 
cufflinks.sh|16|5|73:26:30|3
trinity.sh|20|5|253:54:48|10
GMAP\_build.sh|1|1|00:31:57| 
GMAP\_run.sh|16|5|12:03:48| 
Porticullis|16|1|23:09:12| 
mikado.sh|32|2|73:52:13|4
braker.sh|32|1|184:24:48|8
PASA\_seqclean.sh|16|1|04:29| 
PASA\_alignAssembly.sh|16|1|127:48:58|5
PASA\_annotLoadandCompare.sh|1|50| | 

## Usefull links

For installing and running details, please refer to the following links:

-STAR
https://github.com/alexdobin/STAR

-StringTie
https://github.com/gpertea/stringtie

-PsiCLASS
https://github.com/splicebox/PsiCLASS

-Cufflinks
https://github.com/cole-trapnell-lab/cufflinks

-Trinity
https://github.com/trinityrnaseq/trinityrnaseq/wiki

-GMAP/GSNAP
http://research-pub.gene.com/gmap/

-Portcullis
https://github.com/maplesond/portcullis

-MIKADO
https://github.com/EI-CoreBioinformatics/mikado

-BRAKER
https://github.com/Gaius-Augustus/BRAKER

-PASA
https://github.com/PASApipeline/PASApipeline/wiki

-DESEQ2
https://github.com/mikelove/DESeq2

-CEMiToll
https://github.com/csbl-usp/CEMiTool

-PlantTFDB
http://planttfdb.gao-lab.org

-KOALA
https://www.kegg.jp/blastkoala/
