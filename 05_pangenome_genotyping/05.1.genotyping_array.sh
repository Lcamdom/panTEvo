#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --time=5-02:30:02
#SBATCH --job-name="wild1_pangenome_genotyping"
#SBATCH --array=1-27  # Adjust the range based on the number of populations

module load seqtk/1.3-tstesdq
module load micromamba/1.4.2-7jjmfkf
eval "$(micromamba shell hook --shell=bash)"
micromamba activate vg

# Number of reads to sample
NUM_READS=78000000

# Read the population information from the file
IFS=' ' read -r -a line <<< $(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" populations.txt)
POP=${line[0]}
READ1=${line[1]}
READ2=${line[2]}

# Create a unique copy of the xg file for this job
GRAPH_XG="../GhirPangenome_withTIPs_newindex.xg"
LOCAL_XG="${SLURM_TMPDIR}/${POP}_index.xg"
cp $GRAPH_XG $LOCAL_XG

echo "Processing $POP"

# Subsample reads using seqtk
echo "Subsampling reads for $POP"
seqtk sample -s100 $READ1 $NUM_READS > ${POP}_subsampled_1.fq
seqtk sample -s100 $READ2 $NUM_READS > ${POP}_subsampled_2.fq

# Map reads to the pangenome graph
echo "Mapping reads for $POP"
vg giraffe -t 24 -f ${POP}_subsampled_1.fq -f ${POP}_subsampled_2.fq -x $LOCAL_XG > ${POP}_alignments.gam

# Pack the alignments
echo "Packing alignments for $POP"
vg pack -x $LOCAL_XG -g ${POP}_alignments.gam -o ${POP}_alignments.pack

# Call variants
echo "Calling variants for $POP"
vg call $LOCAL_XG -k ${POP}_alignments.pack > ${POP}.vcf

# Remove the local copy of the xg file
rm $LOCAL_XG

echo "$POP VCF generated"
