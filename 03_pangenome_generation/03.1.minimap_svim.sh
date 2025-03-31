#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --job-name="query_to_REFmapping-SVIM"

## load modules
module load conda

#activate conda environment
conda activate SV-calling

# 1) Align genomes with minimap2
minimap2 -ax asm5 -t 12 AD1_TM1_CRI_chrs.fa \
AD1_CRI12_chrs.fasta > AD1-CRI12_to_AD1-TM1.sam

# 2) SAMtoBAM, sort and index

samtools view -@ 12 -b AD1-CRI12_to_AD1-TM1.sam -o AD1-CRI12_to_AD1-TM1.bam
samtools sort -@ 12 -o AD1-CRI12_to_AD1-TM1_sorted.bam AD1-CRI12_to_AD1-TM1.bam
samtools index AD1-CRI12_to_AD1-TM1_sorted.bam

# 3) detect SVs with SVIM-asm

svim-asm haploid \
--min_sv_size 40 \
./ \
./AD1-CRI12_to_AD1-TM1_sorted.bam \
./AD1_TM1_CRI_chrs.fa

# 4) extract SVs (insertions and deletions)  from vcf file
mv variants.vcf variants_AD1-CRI12_to_AD1-TM1.vcf
grep -E '##|SVTYPE=DEL|SVTYPE=INS' variants_AD1-CRI12_to_AD1-TM1.vcf \
> AD1-CRI12_to_AD1-TM1_indels.vcf
