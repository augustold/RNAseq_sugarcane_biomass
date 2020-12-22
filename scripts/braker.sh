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
