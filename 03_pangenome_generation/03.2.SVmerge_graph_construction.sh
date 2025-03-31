#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=64G
#SBATCH --job-name="build_graph"
#SBATCH --partition=all

# This script merges VCF files and creates a pangenome with VG inside a folder called constructed_

##load environment and modules (samtools needed to run bcftools)
module load conda
module load samtools
module load bcftools

export TMPDIR=./tmp

### first, need to include the following line in the vcf header- make sure it's tab separated:
#CHROM     POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample


## then compress and index vcf files
#ls *indels.vcf > vcffiles.ls
#while read f ; do bgzip "$f" ; done < vcffiles.ls
#ls *vcf.gz > listofvcffiles.txt
#while read f ; do tabix -p vcf "$f" ; done < listofvcffiles.txt

# To start, we merge multiple VCFs (each with their own sample) and ensure there are no multi-allelic entries via:
# Parameters:
#   -m, --merge <string>               allow multiallelic records for <snps|indels|both|all|none|id>#·
#   Error: Duplicate sample names (Sample), use --force-samples to proceed anyway. ---> this is why I added --force-samples
#bcftools merge -m none --force-samples -l listofvcffiles.txt -o bcftoolsmerged.vcf.gz -O z --threads 8

#conda deactivate bcftools
conda activate graph_making

# index bcftools merged
#tabix -p vcf bcftoolsmerged.vcf.gz

# merge highly similar variants with truvari
# truvari collapse to merge similar variants
# -p PCTSEQ, --pctseq PCTSEQ
#                         Min percent sequence similarity. Set to 0 to ignore. (0.95)
# -P PCTSIZE, --pctsize PCTSIZE
#                         Min pct allele size similarity (minvarsize/maxvarsize) (0.95)
#~/.local/bin/./truvari collapse -i bcftoolsmerged.vcf.gz -p 0 -P 0.5 -s 0 -o tru-merged.vcf.gz -c tru-collapsed.vcf

# compress and index the file
#gunzip variants_annotated.vcf.gz
#bgzip variants_annotated.vcf
#tabix -p vcf variants_annotated.vcf.gz

# pangenome construction

# # construct graph wih vg
#vg construct -r ../../B713_to_CRI/AD1_TM1_CRI_chrs.fa \
#-v variants_annotated.vcf.gz \
#-a -S -f -t 8 > GhirPangenome_withTIPs.vg

#vg construct -r /scratch/074-arabidopsis-MITEs/Lana/rice_tips/genomes/Npb.fasta -v tru-merged_rmdups.vcf.gz -a -S -f -t 8 > ricepangenome75_nodupstry.vg


# # store the graph in the xg/gcsa index pair
date
echo '---- Indexing full graph -----'
vg index -t 8 -x GhirPangenome_withTIPs_newindex.xg -k 16 GhirPangenome_withTIPs.vg
echo 'Finished .xg index'
date
vg index -t 8 -Z 4096 -g GhirPangenome_withTIPs_newindex.gcsa -k 16 GhirPangenome_withTIPs.vg
echo 'Finished .gcsa index'

# # subsample graph
#echo subsampling graph
#vg find -x GhirPangenome_newindex.xg -p A01:100000-100500 -E > A01_100000_100500.vg
#vg index -t 8 -x A01_100000_100500.xg -g A01_100000_100500.gcsa -k 16 A01_100000_100500.vg

# # # simplify graph (join 32bp nodes)
#echo simplifying graph
#vg mod -t 8 -u A01_100000_100500.vg -X 256 > A01_100000_100500_mod.vg
#vg view A01_100000_100500_mod.vg > A01_100000_100500_mod.gfa

# # echo indexing simplified graph
#vg index -t 8 -x A01_100000_100500_mod.xg -g A01_100000_100500_mod.gcsa -k 16 A01_100000_100500_mod.vg

#date
echo 'Graph done'
