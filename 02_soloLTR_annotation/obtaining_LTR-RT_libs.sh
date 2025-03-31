#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --job-name="LTR-RT_clustering_by_lineage"

## slurm script to obtain species-level intact LTR-RT libraries with the name convention: "all_${SPECIES}_genomes_intact_LTR_${LINEAGE}_.fasta
### Input files should be all LTR-RTs of a specific lineage for all species

##00. Define your species
SPECIES="my_species_name"

#load dependencies
module load conda
conda activate cd-hit

##01. Sequence clustering at the 80% length and 80% homology
while read -r f; do
    cd-hit-est -i "all_${SPECIES}_genomes_intact_LTR_${f}_.fasta" \
    -o "all_${SPECIES}_genomes_intact_LTR_${f}_cons.fa" \
    -c 0.8 -s 0.8 -T 12 -d 0 -M 48000000
done < "lineages_list"

##02. Concatenating all representatives for every lineage into a single library
cat all_${SPECIES}_genomes*_cons.fa > ${SPECIES}_LTR-RT_lib.fa
