#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=18G
#SBATCH --job-name="getting_lLTRsINTs_libs"

#### Obtaining the lLTR and INT fasta files for each intact LTR-RT from our libraries

genome="my_genome"
species="my_species"

##01. making a table with the coordenates/seqnames of all RT elements
##remove headers and set the order of the fields

grep -v '^#' ${genome}.mod.retriever.all.scn | awk '{printf("%s:%s..%s#LTR/%s\t%s-%s\t%s-%s\n", $12, $1, $2, $19, 1, $6+1, $6+2, $7-$1 )}' > ${genome}_lLTR_INTs_coords.txt

##do in a loop for all files (only if you're doing it for many files at once, where all my accession names are in a text file "${species}_accessions"):

#while read f ; do grep -v '^#' ${species}_"$f"_chrs.fa.mod.retriever.all.scn | 
#awk '{printf("%s:%s..%s#LTR/%s\t%s-%s\t%s-%s\n", $12, $1, $2, $19, 1, $6+1, $6+2, $7-$1 )}' \
#> ${species}_"$f"_chrs.fa.mod.retriever_noheader.scn ; done < ${species}_accessions
#cat ${species}*_noheader.scn > ${species}_genomes_lLTR_INTs_coords.txt

#get only coordinates to grep the reps from the other file

grep ">" ${genome}_LTR-RTlib.fa | sed 's/>//' > ${genome}_LTR-RTlib_names.ls

sed 's/.*:\([^#]*\)#.*/\1/g' ${genome}_LTR-RTlib_names.ls > ${genome}_LTR-RTlib_names_mod.ls

while read f ; do grep "$f" ${genome}_lLTR_INTs_coords.txt ; done < ${genome}_LTR-RTlib_names_mod.ls > ${genome}_lLTR_INTs_coords_cons.txt
sort -u ${genome}_lLTR_INTs_coords_cons.txt > ${genome}_lLTR_INTs_coords_cons_SU.txt

### making a two-column file to get the full headers and the "match" strings next to each other

sed 's/.*://' ${genome}_LTR-RTlib_names.ls > ${genome}_LTR-RTlib_names_1.ls
sed 's/_[^_]*$//' ${genome}_LTR-RTlib_names_1.ls > ${genome}_LTR-RTlib_names_1_2.ls
paste ${genome}_LTR-RTlib_names.ls ${genome}_LTR-RTlib_names_1_2.ls > ${genome}_LTR-RTlib_names_lib_2columns.ls

sed 's/.*://' ${genome}_lLTR_INTs_coords_cons_SU.txt > ${genome}_lLTR_INTs_coords_cons_SU_mod.txt

## find correspondances between files

awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1]=$2"\t"$3; next} {print $0, a[$2]}' ${genome}_lLTR_INTs_coords_cons_SU_mod.txt \
 ${genome}_LTR-RTlib_names_lib_2columns.ls > ${genome}_intact_names_lib_with_lLTR_INT_coords.txt


## print "region" files for lLTRs and INTs

awk '{printf("%s:%s\n", $1, $3)}' ${genome}_intact_names_lib_with_lLTR_INT_coords.txt > ${genome}_intact_lib_lLTR_regfile.txt

awk '{printf("%s:%s\n", $1, $4)}' ${genome}_intact_names_lib_with_lLTR_INT_coords.txt > ${genome}_intact_lib_INTs_regfile.txt

## run samtools to get the fastas

module load samtools
samtools faidx ${species}_LTR-RTlib.fa
samtools faidx -r ${genome}_intact_lib_INTs_regfile.txt G.raimondii_LTR-RTlib.fa > ${genome}_intact_LTRlib.fa
samtools faidx -r ${genome}_intact_lib_lLTR_regfile.txt G.raimondii_LTR-RTlib.fa > ${genome}_intact_INTlib.fa
