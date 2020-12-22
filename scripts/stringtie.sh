cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load StringTie/1.3.5

date

stringtie /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/atlas_merged.bam \
--rf -p 16 -o /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/stringtie/atlas_ST.gtf

date

