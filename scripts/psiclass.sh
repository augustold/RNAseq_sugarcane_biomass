cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/psiclass

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1

date

/sonas-hs/it/hpc/home/bnb/src/psiclass/psiclass/psiclass -b /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/atlas_merged.bam -p 16 -o /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/psiclass

mv psiclass_* psiclass/

date
