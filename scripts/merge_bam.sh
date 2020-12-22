cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load SAMtools/1.9

date

samtools merge atlas_merged.bam aligned_2/*.bam
samtools index atlas_merged.bam

date

