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
