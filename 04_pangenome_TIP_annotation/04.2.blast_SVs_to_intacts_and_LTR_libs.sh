#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --export=NONE
#SBATCH --job-name="blast_SVs"


module load blast

makeblastdb -in Ghir_tru-merged_INDELs_filtered.fasta \
-input_type fasta \
-dbtype nucl

blastn -query all_AD1_genomes_intact_lib.fa \
-db Ghir_tru-merged_INDELs_filtered.fasta \
-outfmt "6 qseqid sseqid qstart qend sstart send pident length qlen slen qcovs" \
-perc_identity 80 \
-qcov_hsp_perc 80 \
-num_threads 12 \
-out Ghir_tru-merged_INDELs_intacts_blasthits.tsv

blastn -query all_AD1_genomes_intact_lib_lLTR.fa \
-db Ghir_tru-merged_INDELs_filtered.fasta \
-outfmt "6 qseqid sseqid qstart qend sstart send pident length qlen slen qcovs" \
-perc_identity 80 \
-qcov_hsp_perc 80 \
-num_threads 12 \
-out Ghir_tru-merged_INDELs_LTRs_blasthits.tsv


exit
