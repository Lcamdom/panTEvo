setwd("/Users/lcampos/Desktop/cotton_data/eQTL_work/RNA-seq")

# Load packages
library(DESeq2)

# Read your counts matrix
counts <- read.table("salmon.isoform.counts.matrix", header = TRUE, row.names = 1, sep = "\t")

# Check data
head(counts)

# Create dummy colData with one condition (e.g., "group") just to satisfy DESeq2
sample_info <- data.frame(row.names = colnames(counts), condition = rep("A", ncol(counts)))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = sample_info, design = ~ 1)

# VST transformation
vst_counts <- varianceStabilizingTransformation(dds, blind = TRUE)

# Extract matrix
vst_matrix <- assay(vst_counts)

# Check result
head(vst_matrix)

write.table(vst_matrix, file = "vst_matrix.tsv", sep = "\t", quote = FALSE)