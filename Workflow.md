Create all folders that will be used in this workflow
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280
mkdir scripts star_index star_index_2ndPass aligned_1 aligned_2 stringtie cufflinks trinity psiclass mikado
```

This is the path for script files
```
/sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/scripts
```

### STAR genome index for 1-pass

script: star_1stpass.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load STAR/2.7.0d

date

STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/star_index \
--genomeFastaFiles /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
--sjdbGTFfile /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg_scga7.sort.gff3 \
--genomeChrBinNbits 13 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 150 \
--limitGenomeGenerateRAM 350000000000

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G star_1stpass.sh
```

### Align reads for 1-pass

script: star_align_1stpass.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load STAR/2.7.0d

date

for i in $(ls /sonas-hs/ware/hpc_norepl/data/diniz/fastq/SP80 | sed s/_[12].fq.gz// | sort -u)
do
    STAR \
    --genomeDir /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/star_index \
    --readFilesIn /sonas-hs/ware/hpc_norepl/data/diniz/fastq/SP80/${i}_1.fq.gz,/sonas-hs/ware/hpc_norepl/data/diniz/fastq/SP80/${i}_2.fq.gz \
    --runThreadN 16 \
    --outFileNamePrefix aligned_1/${i}. \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes All \
    --chimSegmentMin 12 \
    --chimSegmentReadGapMax 3 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterIntronMotifs RemoveNoncanonical \
    --clip5pNbases 0 \
    --seedSearchStartLmax 50 \
    --genomeLoad LoadAndKeep \
    --limitBAMsortRAM 35000000000
done

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G star_align_1stpass.sh
```

### STAR genome index for 2-pass

script: star_2ndpass.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/aligned_1/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load STAR/2.7.0d

date

STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/star_index_2ndPass \
--genomeFastaFiles /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
--sjdbGTFfile /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg_scga7.sort.gff3 \
--genomeChrBinNbits 13 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 150 \
--limitGenomeGenerateRAM 350000000000 \
--limitSjdbInsertNsj 2000000 \
--sjdbFileChrStartEnd SP80L111.SJ.out.tab SP80L121.SJ.out.tab SP80L131.SJ.out.tab SPI11.SJ.out.tab SPI12.SJ.out.tab SPI13.SJ.out.tab SPI51.SJ.out.tab SPI52.SJ.out.tab SPI53.SJ.out.tab SPI91.SJ.out.tab SPI92.SJ.out.tab SPI93.SJ.out.tab SPL11.SJ.out.tab SPL12.SJ.out.tab SPL13.SJ.out.tab SPR1.SJ.out.tab SPR2.SJ.out.tab SPR3.SJ.out.tab

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G star_2ndpass.sh
```

### Align reads for 2nd-pass

script: star_align_2ndpass.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load STAR/2.7.0d

date

for i in $(ls /sonas-hs/ware/hpc_norepl/data/diniz/fastq/SP80 | sed s/_[12].fq.gz// | sort -u)
do
    STAR \
    --genomeDir /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/star_index_2ndPass \
    --readFilesIn /sonas-hs/ware/hpc_norepl/data/diniz/fastq/SP80/${i}_1.fq.gz,/sonas-hs/ware/hpc_norepl/data/diniz/fastq/SP80/${i}_2.fq.gz \
    --runThreadN 16 \
    --outFileNamePrefix aligned_2/${i}. \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes All \
    --chimSegmentMin 12 \
    --chimSegmentReadGapMax 3 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterIntronMotifs RemoveNoncanonical \
    --clip5pNbases 0 \
    --seedSearchStartLmax 50 \
    --genomeLoad LoadAndKeep \
    --limitBAMsortRAM 35000000000
done

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G star_align_2ndpass.sh
```

### Merge BAM files

script: merge_bam.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load SAMtools/1.9

date

samtools merge atlas_merged.bam aligned_2/*.bam
samtools index atlas_merged.bam

date
```
```
qsub -cwd -pe threads 1 -l m_mem_free=20G merge_bam.sh
```

### Genome-guided: StringTie

script: stringtie.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/stringtie

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load StringTie/1.3.5

date

stringtie ../atlas_merged.bam \
--rf -p 16 -o atlas_ST.gtf

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G stringtie.sh
```

### Genome-guided: psyCLASS

script: psiclass.sh

```
#!/bin/bash

cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1

date

/sonas-hs/it/hpc/home/bnb/src/psiclass/psiclass/psiclass -b /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/atlas_merged.bam -p 16 -o /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/psiclass

mv psiclass_* psiclass/

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G psiclass.sh
```

### Genome-guided: Cufflinks

script: cufflinks.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/

date

/sonas-hs/ware/hpc/home/bwang/software/cufflinks-2.1.1.Linux_x86_64/cufflinks \
-p 16 \
-o /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/cufflinks \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg_scga7.sort.gff3 \
--library-type fr-firststrand \
-u /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/atlas_merged.bam

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G cufflinks.sh
```

### Genome-guided: Trinity

script: trinity.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load SAMtools/1.9 
module load Jellyfish/2.2.10
module load Salmon/0.14.1
module load Trinity/2.8.4

date

Trinity \
--genome_guided_bam /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/atlas_merged.bam \
--genome_guided_max_intron 10000 \
--max_memory 100G --CPU 20 \
--output /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/trinity/trinity_out_dir

date
```
```
qsub -cwd -pe threads 20 -l m_mem_free=5G trinity.sh
```

### Align Trinity transcripts back to the genome using GMAP-GSNAP

#### Build gmap index

script: GMAP_build.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load GMAP-GSNAP/2019-03-15
module load SAMtools/1.9

date

gmap_build -d gmap_index -D . -k 13 /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa

date
```
```
qsub -cwd -pe threads 1 -l m_mem_free=16G GMAP_build.sh
```

#### Alignment

script: GMAP_run.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/trinity/trinity_out_dir

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load GMAP-GSNAP/2019-03-15
module load SAMtools/1.9

date

gmap \
-D /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280 \
-d gmap_index \
-f gff3_gene \
-t 16 \
-n 1 \
--min-trimmed-coverage=0.70 \
--min-identity=0.95 \
Trinity-GG.fasta > trinity.gff3

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G GMAP_run.sh
```

### PORTCULLIS

Atention: run this step on brie.cshl.edu

Copy bam and genome reference files from bnbdev.cshl.edu to brie.cshl.edu
```
scp -r /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/atlas_merged.bam diniz@brie.cshl.edu:/home/diniz/SP80-3280/
scp -r /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa diniz@brie.cshl.edu:/home/diniz/SP80-3280/
```

Run Portcullis
```
screen
cd /home/diniz/SP80-3280
#Before run Portcullis, activate conda enviroment where it is installed
source activate portcullis
 
date

portcullis full \
   --threads 16 \
   --output portcullis_out \
   sc.mlc.cns.sgl.utg.scga7.importdb.fa atlas_merged.bam
   
date
```

### MIKADO

Prepare 'mikado' directory

```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/mikado

# copy portcullis junction.bed file
cp /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/portcullis_out/3-filt/portcullis_filtered.pass.junctions.bed .

# get uniprot_sprot_plants.fasta
cp /sonas-hs/ware/hpc_norepl/data/xwang/NAM/Sorg_newWF/mikado_class2/uniprot_sprot_plants.fasta .

# copy an rename gtf/gff to directory

# cufflinks
cp /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/cufflinks/transcripts.gtf .
mv transcripts.gtf cufflinks.gtf

# stringtie
cp /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/stringtie/atlas_ST.gtf .
mv atlas_ST.gtf stringtie.gtf

# stringtie
cp /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/stringtie/atlas_ST.gtf .
mv atlas_ST.gtf stringtie.gtf

# trinity
cp /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/trinity/trinity_out_dir/trinity.gff3 .
```

Prepare list.txt file

```
cufflinks.gtf	cuff	True		False
psiclass.gtf	cl	False		False
stringtie.gtf	st	True	1	False
trinity.gff3	tr	False	-0.5	False
```

### Running MIKADO pipeline

script: mikado.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/mikado

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6  
module load TransDecoder/5.5.0-Perl-5.28.0

date

#Creating the configuration file for Mikado
mikado configure \
--list list.txt \
--reference /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-t 32 \
--mode permissive \
--scoring plant.yaml  \
--copy-scoring plant.yaml \
--junctions portcullis_filtered.pass.junctions.bed \
-bt uniprot_sprot_plants.fasta \
configuration.yaml

#Mikado prepare
mikado prepare --json-conf configuration.yaml

#BLAST of the candidate transcripts
makeblastdb -in uniprot_sprot_plants.fasta -dbtype prot -parse_seqids > blast_prepare.log

blastx -max_target_seqs 5 -num_threads 32 -query mikado_prepared.fasta -outfmt 5 -db uniprot_sprot_plants.fasta -evalue 0.000001 2> blast.log | sed '/^$/d' | gzip -c - > mikado.blast.xml.gz

TransDecoder.LongOrfs -t mikado_prepared.fasta
TransDecoder.Predict -t mikado_prepared.fasta

#Mikado serialise
mikado serialise --json-conf configuration.yaml --xml mikado.blast.xml.gz --orfs mikado_prepared.fasta.transdecoder.bed --blast_targets uniprot_sprot_plants.fasta --transcripts mikado_prepared.fasta

#Mikado pick
mikado pick --json-conf configuration.yaml --subloci-out mikado.subloci.gff3 --procs 32

gffread mikado.loci.gff3 \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-x mikado.loci.cds.fasta -y mikado.loci.protein.fasta

date
```
```
qsub -cwd -pe threads 32 -l m_mem_free=2G mikado.sh
```

# Insert Mikado downstream analysis here

### BRAKER

Using 'Hints' instead of 'bam'

Convert the bam files to RNA-seq hints file 
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load SAMtools/1.9

date

samtools sort atlas_merged.bam -o atlas_merged_sorted.bam

/sonas-hs/ware/hpc/home/diniz/software/Augustus/bin/bam2hints --intronsonly --in=atlas_merged_sorted.bam --out=SP80-3280.RNAseq.hints

date
```

script: braker.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280

module load foss/2018b
module load Python/3.6.6

date

/sonas-hs/ware/hpc/home/diniz/software/BRAKER/scripts/braker.pl \
--cores 32 \
--min_contig=10000 \
--DIAMOND_PATH=/sonas-hs/ware/hpc/home/diniz/software/diamond \
--genome=/sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
--hints=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/SP80-3280.RNAseq.hints

date
```
```
qsub -cwd -pe threads 32 -l m_mem_free=1G braker.sh
```

script: braker_masked.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280

module load foss/2018b
module load Python/3.6.6

date

/sonas-hs/ware/hpc/home/diniz/software/BRAKER/scripts/braker.pl \
--cores 32 \
--min_contig=10000 \
--DIAMOND_PATH=/sonas-hs/ware/hpc/home/diniz/software/diamond \
--genome=/sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
--hints=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/SP80-3280.RNAseq.hints

date
```
```
qsub -cwd -pe threads 32 -l m_mem_free=1G braker_masked.sh
```

### BRAKER - RNAseq hints + mikado proteins

Prepare mikado proteins fasta

```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/mikado_3rd_test

#Simplify fasta headers
sed -e 's/ .*//g' mikado.loci.protein.fasta > mikado.loci.protein_newheader.v1.fasta
sed -r 's/\s+//g' mikado.loci.protein_newheader.v1.fasta > mikado.loci.protein_newheader.v2.fasta
sed -r 's/\.//g' mikado.loci.protein_newheader.v2.fasta > mikado.loci.protein_newheader.v3.fasta

#Run Pfam fo filter proteins
#Get proteins with domain

/sonas-hs/ware/hpc_norepl/data/programs/maker_v3_updated/maker/bin/fasta_tool --select mikado.protein.pfam.list.txt mikado.loci.protein_newheader.v3.fasta > mikado.loci.protein_newheader.v3.Pfam.fasta

```

### ATENTION!

```diff
! WARNING: in file /sonas-hs/ware/hpc/home/diniz/software/BRAKER/scripts/braker.pl at line 1256 file /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker/genome.fa contains a highly fragmented assembly (450608 scaffolds). This may lead to problems when running AUGUSTUS via braker in parallelized mode. You set --cores=16. You should run braker.pl in linear mode on such genomes, though (--cores=1).
```

script: braker_prot.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280

module load foss/2018b
module load Python/3.6.6

date

/sonas-hs/ware/hpc/home/diniz/software/BRAKER/scripts/braker.pl \
--cores 1 \
--min_contig=10000 \
--DIAMOND_PATH=/sonas-hs/ware/hpc/home/diniz/software/diamond \
--genome=/sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
--hints=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/SP80-3280.RNAseq.hints \
--AUGUSTUS_ab_initio \
--prot_seq=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/mikado_3rd_test/mikado.loci.protein_newheader.v4.Pfam.fasta \
--prg=gth --gth2traingenes

date
```
```
qsub -cwd -pe threads 1 -l m_mem_free=32G braker_prot.sh
```

# PASA tool for annotation update.

## Install PASA

```
# Prepare working directory

cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis
mkdir PASA_run; cd PASA_run; mkdir scripts

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
module load git/2.19.1

conda create -n pasa
source activate pasa
conda install -c bioconda pasa
conda install -c bioconda/label/cf201901 gmap
conda install -c bioconda/label/cf201901 blat
conda install openssl=1.0
conda install -c bioconda/label/cf201901 perl-dbi
conda install -c bioconda perl-dbd-sqlite
conda install -c bioconda perl-uri
conda install -c bioconda perl-db-file

wget https://github.com/PASApipeline/PASApipeline/releases/download/pasa-v2.4.1/PASApipeline.v2.4.1.FULL.tar.gz

tar -zxvf PASApipeline.v2.4.1.FULL.tar.gz
rm -rf PASApipeline.v2.4.1.FULL.tar.gz
cd PASApipeline.v2.4.1; make
```

## Input
- Repeatmasked genome
- Trasncripts
	- ESTs
	- Full-Lengths (Nishiyama 2014)
	- Mikado CDS

## Step 1: combine transcripts

```
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280

cat \
ESTs_SP80-3280.fasta \
sugarcane.fulllength.analysis.all.fasta \
/sonas-hs/ware/hpc_norepl/data/diniz/analysis/mikado_3rd_test/mikado.loci.TErmv.cds.fasta \
> /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run/SP80.est.flc.mikado.combined.fasta
```

Also get the full-lengths IDs in a different file
```
grep -e ">" sugarcane.fulllength.analysis.all.fasta | sed 's/>//' > /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run/FL_accs.txt
```

## Step 2: cleaning the transcript sequences

script: seqclean.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run
 
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate pasa

date

PASApipeline.v2.4.1/bin/seqclean SP80.est.flc.mikado.combined.fasta -c 16

mkdir cleaning
mv cleaning_* cleaning/

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=1G seqclean.sh 
```
  
Output summary:

```
less seqcl_SP80.est.flc.mikado.combined.fasta.log

**************************************************
Sequences analyzed:    577771
-----------------------------------
                   valid:    577694  (5503 trimmed)
                 trashed:        77
**************************************************
----= Trashing summary =------
            by 'low_qual':        4
               by 'short':       33
              by 'shortq':       24
                by 'dust':       16
------------------------------
Output file containing only valid and trimmed sequences: SP80.est.flc.mikado.combined.fasta.clean
For trimming and trashing details see cleaning report  : SP80.est.flc.mikado.combined.fasta.cln
--------------------------------------------------
```

## Step 3: Transcript alignments followed by alignment assembly

Before running the transcript alignment step, copy and edit the copy files annotCompare.config and alignAssembly.config from PASA installation folder.

```
#cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run

cp /sonas-hs/ware/hpc_norepl/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/B97/alignAssembly.config .
cp /sonas-hs/ware/hpc_norepl/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/B97/annotCompare.config .
```

Edit the database path :

```
# database settings
DATABASE=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run/sqlite/SP80.sqlite
```

Create a sqlite folder and set the permission

```
mkdir sqlite; chmod -R 777 sqlite
```

Run the PASA alignment assembly pipeline

script: alignAssembly.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run
 
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
module load GMAP-GSNAP/2019-03-15

source activate pasa

date

PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c alignAssembly.config \
-C -R \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
-t SP80.est.flc.mikado.combined.fasta.clean -T \
-u SP80.est.flc.mikado.combined.fasta \
-f FL_accs.txt \
--ALIGNERS gmap \
--CPU 16

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh 
```

## Step 4: Annotation Comparisons and Annotation Updates

Loading your preexisting protein-coding gene annotations and performing an annotation comparison and generating an updated gene set

script: annotLoadandCompare.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run 
 
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate pasa

date

PASApipeline.v2.4.1/scripts/Load_Current_Gene_Annotations.dbi \
-c alignAssembly.config \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
-P /sonas-hs/ware/hpc_norepl/data/diniz/analysis/mikado_3rd_test/mikado.loci.TErmv.gff3

PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c annotCompare.config \
-A \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
--CPU 16 -t SP80.est.flc.mikado.combined.fasta.clean

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=3G annotLoadandCompare.sh 
```

## Step 5: consider run this
```
https://github.com/PASApipeline/PASApipeline/wiki/PASA_abinitio_training_sets
```

script: abinitiotraining.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run 
 
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate pasa

date

PASApipeline.v2.4.1/scripts/pasa_asmbls_to_training_set.dbi \
--pasa_transcripts_fasta SP80.sqlite.assemblies.fasta \
--pasa_transcripts_gff3 SP80.sqlite.pasa_assemblies.gff3

date
```
```
qsub -cwd -pe threads 1 -l m_mem_free=16G abinitiotraining.sh
```
