#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --job-name="find_and_classify_intactLTR-RTs"

#### Slurm script to be run for all genome assemblies to find intact LTR-RTs and classify these into lineages using TEsorter

#0. Define input files and output prefix
GENOME_FASTA="my_genome_chrs.fa"
INTACTS_FASTA="my_genome_chrs.fa.mod.LTR.intact.fa"
OUTPUT_PREFIX="my_genome_intacts_TEsorter"


#1. EDTA raw LTR module to find intact LTR-RTs

##loading dependencies
module load conda
conda activate EDTA

##EDTA command
EDTA_raw.pl --genome $GENOME_FASTA --type ltr

deactivate

#2. TEsorter on the intacts.fa EDTA output to further classify them into lineages

##loading dependencies
module load conda
conda activate TEsorter

##TEsorter command
TEsorter  $INTACTS_FASTA \
-db rexdb-plant \
-pre $OUTPUT_PREFIX
