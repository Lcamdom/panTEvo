#!/bin/bash
#SBATCH --mem-per-cpu=18G
#SBATCH --job-name="TIPannotation"

##first: merge soloLTR and intactLTR-RT bedfiles

#specify solo and intact LTR tips:

#sed 's/#LTR/#soloLTR/g' Ghir_tru-merged_INDELs_LTRs_blasthits.tsv_0.8SV_corrected_sorted_merged.bed > Ghir_tru-merged_INDELs_LTRs_blasthits.tsv_0.8SV_corrected_sorted_merged_mod.bed

#sed 's/#LTR/#LTR-RT/g'Ghir_tru-merged_INDELs_intacts_blasthits.tsv_0.8SV_corrected_sorted_merged.bed > Ghir_tru-merged_INDELs_intacts_blasthits.tsv_0.8SV_corrected_sorted_merged_mod.bed

#cat Ghir_tru-merged_INDELs_intacts_blasthits.tsv_0.8SV_corrected_sorted_merged_mod.bed \
#Ghir_tru-merged_INDELs_LTRs_blasthits.tsv_0.8SV_corrected_sorted_merged_mod.bed > Ghir_tru-merged_all_TIPs.bed

module load bedtools
#
#bedtools sort -i Ghir_tru-merged_all_TIPs.bed \
#> Ghir_tru-merged_all_TIPs_sorted.bed


# Annotating the VCF file

# 1. Obtaining variant ID + TE annotation file
awk '{print $1,$5}' Ghir_tru-merged_all_TIPs_sorted.bed > temp_id_info.txt

awk '{print $1,";TE_ANNOT="$2}' temp_id_info.txt > id_info.txt

head id_info.txt

# 2. Add the INFO line to the header
sed '/##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">/a ##INFO=<ID=TE_ANNOT,Number=1,Type=String,Description="Transposable element annotation">' tru-merged.vcf > temp_variants.vcf

# 3. Separate the header
grep '#' temp_variants.vcf > temp_header.txt

# 4. Separate the variant lines
grep -v "#" temp_variants.vcf > temp_variants_nohead.vcf

# 5. Add TE annotation
#awk 'FNR==NR{a[$1]=$2;next}{if(a[$3]==""){a[$3]=na}; printf "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n",$1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7,"\t",$8,a[$3],"\t",$9,"\t",$10, "\t",$11,"\t",$12,"\t",$13,"\t",$14,"\t",$15,"\t",$16}' \
#id_info.txt temp_variants_nohead.vcf > temp_variants_annot.vcf

awk '
FNR == NR {
    a[$1] = $2;
    next
}
{
    # Split the current line into fields using tab as the delimiter
    split($0, fields, "\t");

    # Check if the annotation exists, if not use "na"
    if (a[fields[3]] == "") {
        annotation = "na";
    } else {
        annotation = a[fields[3]];
    }

    # Print the fields with the annotation inserted after the 8th field
    for (i = 1; i <= 8; i++) {
        printf "%s\t", fields[i];
    }
    printf "%s\t", annotation;
    for (i = 9; i <= length(fields); i++) {
        printf "%s%s", fields[i], (i == length(fields) ? "\n" : "\t");
    }
}' id_info.txt temp_variants_nohead.vcf > temp_variants_annot.vcf


# 6. Concatenate header and variants
cat temp_header.txt temp_variants_annot.vcf > variants_annotated.vcf

