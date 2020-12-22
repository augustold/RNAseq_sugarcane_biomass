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
