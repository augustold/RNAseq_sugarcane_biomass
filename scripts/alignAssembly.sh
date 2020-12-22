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
