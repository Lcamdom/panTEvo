#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=52G
#SBATCH --job-name="soloLTRs"

##### Slurm script to find soloLTRs in a genome. It will mask the genome for intact LTR-RTs, blast the LTR and INT libraries to the genome, filter the results, and only keep the LTRhits with no INT regions 1kbp around. 

### 0. name your input 
accession="my_accession"

#### 1. masking genome fastas
module load bedtools

cat ${accession}_chrs.fa.mod.harvest.combine.gff3 \
${accession}_chrs.fasta.mod.LTR.intact.gff3 | \
bedtools sort | bedtools merge -i - \
> ${accession}.chrs.fasta.intacts.combined.gff3

bedtools maskfasta \
-fi ../${accession}_chrs.fa \
-bed ${accession}.chrs.fa.intacts.combined.gff3 \
-fo ${accession}_chrs_maskedRT-LTRs.fa


### 2. BLASTing LTRs and INTs

module load blast

makeblastdb -in ${accession}_chrs_maskedRT-LTRs.fa \
-input_type fasta \
-dbtype nucl

blastn -query ${accession}_intact_LTRlib.fa \
-db ${accession}_chrs_maskedRT-LTRs.fa \
-outfmt "6 qseqid sseqid qstart qend sstart send pident length qlen slen qcovs" \
-perc_identity 80 \
-qcov_hsp_perc 80 \
-num_threads 8 \
-out ${accession}_chrs_vs_LTRlib_blasthits1.tsv

blastn -query ${accession}_intact_INTlib.fa \
-db ${accession}_chrs_maskedRT-LTRs.fa \
-outfmt "6 qseqid sseqid qstart qend sstart send pident length qlen slen qcovs" \
-perc_identity 80 \
-qcov_hsp_perc 20 \
-num_threads 8 \
-out ${accession}_chrs_vs_INTlib_blasthits.tsv


###3. finding_solos

## 3.1 modify the outputs for bedtools:
##3.1.1 LTRs

awk 'BEGIN{OFS="\t"} \
{print $2, $5-1, $6-1, $1, "."}' \
${accession}_chrs_vs_LTRlib_blasthits.tsv |\
awk -v OFS='\t' '{ if ($2 > $3) \
{ temp = $2; $2 = $3; $3 = temp; orientation = "-" } \
else { orientation = "+" } print $0, orientation }' \
> ${accession}_chrs_vs_LTRlib_blasthits_corrected.bed

bedtools sort -i ${accession}_chrs_vs_LTRlib_blasthits_corrected.bed \
> ${accession}_chrs_vs_LTRlib_blasthits_corrected_sorted.bed

bedtools merge -c 4,4 -o count,first \
-i ${accession}_chrs_vs_LTRlib_blasthits_corrected_sorted.bed \
> ${accession}_chrs_vs_LTRlib_blasthits_corrected_sorted_merged.bed

##3.1.2 INTs

awk 'BEGIN{OFS="\t"} \
{print $2, $5-1, $6-1, $1, "."}' \
${accession}_chrs_vs_INTlib_blasthits.tsv |\
awk -v OFS='\t' '{ if ($2 > $3) \
{ temp = $2; $2 = $3; $3 = temp; orientation = "-" } \
else { orientation = "+" } print $0, orientation }' \
> ${accession}_chrs_vs_INTlib_blasthits_corrected.bed

bedtools sort -i ${accession}_chrs_vs_INTlib_blasthits_corrected.bed \
> ${accession}_chrs_vs_INTlib_blasthits_corrected_sorted.bed

bedtools merge -c 4,4 -o count,first \
-i ${accession}_chrs_vs_INTlib_blasthits_corrected_sorted.bed \
> ${accession}_chrs_vs_INTlib_blasthits_corrected_sorted_merged.bed


## 3.2 get upstream and downstream regs

awk -v OFS='\t' \
'{ print $1, $2 - 1000, $2, $2, $3, $5 }' \
${accession}_chrs_vs_LTRlib_blasthits_corrected_sorted_merged.bed \
> ${accession}_chrs_vs_LTRlib_blasthits_1kbupstream.bed

awk -v OFS='\t' \
'{ print $1, $3, $3 + 1000, $2, $3, $5 }' \
${accession}_chrs_vs_LTRlib_blasthits_corrected_sorted_merged.bed \
> ${accession}_chrs_vs_LTRlib_blasthits_1kbdownstream.bed

## 3.3 filter those out that overlap with INTs

bedtools subtract -A \
-a ${accession}_chrs_vs_LTRlib_blasthits_1kbupstream.bed \
-b ${accession}_chrs_vs_INTlib_blasthits_corrected_sorted_merged.bed \
> ${accession}_chrs_vs_LTRlib_blasthits_1kbupstream_noINTs.bed

 bedtools subtract -A \
 -a ${accession}_chrs_vs_LTRlib_blasthits_1kbdownstream.bed \
 -b ${accession}_chrs_vs_INTlib_blasthits_corrected_sorted_merged.bed \
 > ${accession}_chrs_vs_LTRlib_blasthits_1kbdownstream_noINTs.bed

 ## 3.4 get the original coordinates and only the LTRs that don't have INTs at either side

 awk -v OFS='\t' '{print $1, $4, $5, $6}' \
 ${accession}_chrs_vs_LTRlib_blasthits_1kbdownstream_noINTs.bed \
 > ${accession}_chrs_vs_LTRlib_blasthits_1kbdownstream_noINTs_originalregs.bed

 awk -v OFS='\t' '{print $1, $4, $5, $6}' \
${accession}_chrs_vs_LTRlib_blasthits_1kbupstream_noINTs.bed \
 > ${accession}_chrs_vs_LTRlib_blasthits_1kbupstream_noINTs_originalregs.bed

bedtools intersect \
-a ${accession}_chrs_vs_LTRlib_blasthits_1kbupstream_noINTs_originalregs.bed \
-b ${accession}_chrs_vs_LTRlib_blasthits_1kbdownstream_noINTs_originalregs.bed \
> ${accession}_chrs_vs_LTRlib_blasthits_noINTs_1kbaround.bed
