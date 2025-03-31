#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=18G
#SBATCH --job-name="filtering_TIPs"

hitfile="Ghir_tru-merged_INDELs_intacts_blasthits.tsv"

# filter out everything where hit ($8) =or > 0.8SV ($10)
awk -F'\t' '$8 > 0.8 * $10 {print $0}' ${hitfile} > ${hitfile}_0.8SV.tsv

## 2. modify the outputs for bedtools:

module load bedtools
awk 'BEGIN{OFS="\t"} \
{print $2, $5-1, $6-1, $1, "."}' \
${hitfile}_0.8SV.tsv |\
awk -v OFS='\t' '{ if ($2 > $3) \
{ temp = $2; $2 = $3; $3 = temp; orientation = "-" } \
else { orientation = "+" } print $0, orientation }' \
> ${hitfile}_0.8SV_corrected.bed

bedtools sort -i ${hitfile}_0.8SV_corrected.bed \
> ${hitfile}_0.8SV_corrected_sorted.bed

bedtools merge -c 4,4 -o count,first \
-i ${hitfile}_0.8SV_corrected_sorted.bed \
> ${hitfile}_0.8SV_corrected_sorted_merged.bed
