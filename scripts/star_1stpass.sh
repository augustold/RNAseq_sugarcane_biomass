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
