#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --export=NONE
#SBATCH --job-name="SV_linearising_and_filtering"

### extract fasta seqs from vcf
###1. Remove lines with '##' (header)
grep -v '##' Ghir_tru-merged.vcf | \
### 2. extract columns 3 and 4 (seq name and seq) for deletions
grep "DEL" |\
cut -f3,4 |\
### 3. Change tabs for new lines to make them fasta format
tr "\t" "\n" |\
### 4. Add ">" to the lines starting with "svim" and print it into a fasta
sed -e '/^[0-9]/ s/./>&/' > Ghir_tru-merged_DELs.fasta

##now extract insertions
grep -v '##' Ghir_tru-merged.vcf |\
grep "INS" | cut -f3,5 | tr "\t" "\n" |\
sed -e '/^[0-9]/ s/./>&/' > Ghir_tru-merged_INS.fasta

##concatenate fastas
cat Ghir_tru-merged_DELs.fasta Ghir_tru-merged_INS.fasta \
> Ghir_tru-merged_INDELs.fasta
#rm Ghir_tru-merged_DELs.fasta
#rm Ghir_tru-merged_INS.fasta

perl /scratch/075-melo-TEmovement/COTTON/lucia_panTEs/test/with_curated_genomes/prinseq/./prinseq-lite.pl \
-fasta Ghir_tru-merged_INDELs.fasta \
-min_len 100

mv Ghir_tru-merged_INDELs_prinseq_good*.fasta Ghir_tru-merged_INDELs_filtered.fasta
