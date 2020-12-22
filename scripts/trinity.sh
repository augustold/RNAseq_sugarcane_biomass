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
