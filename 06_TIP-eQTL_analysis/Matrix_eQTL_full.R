# If not installed
install.packages("MatrixEQTL")

# Load it
library(MatrixEQTL)
setwd("/Users/lcampos/Desktop/cotton_data/eQTL_work/50_var_16DPM")


# Load packages
library(MatrixEQTL)

genotype_file <- "transformed_genotype_matrix_TIPsonly.tsv"
expression_file <- "vst_matrix.tsv"

# === Step 1: Load raw data as data.frames for filtering ===

# Load genotype and expression data as matrices first (before converting to SlicedData)
geno_df <- read.table(genotype_file, header = TRUE, row.names = 1, check.names = FALSE)
expr_df <- read.table(expression_file, header = TRUE, row.names = 1, check.names = FALSE)

# MAF Filter Function

# MAF Filter Function - 1st Filter
sv_maf_filter1 <- function(geno_row, threshold = 0.05) {
  # Remove missing values
  geno_row <- geno_row[!is.na(geno_row)]
  
  # Check if the row still has data after removing NAs
  if (length(geno_row) == 0) return(FALSE)
  
  # Count all alleles (0, 1, 2) — approximate allele counts
  counts <- table(geno_row)
  
  # Ensure counts for each genotype (0, 1, 2) are defined
  allele0_count <- ifelse("0" %in% names(counts), counts["0"] * 2, 0) +
    ifelse("1" %in% names(counts), counts["1"] * 1, 0)
  allele1_count <- ifelse("2" %in% names(counts), counts["2"] * 2, 0) +
    ifelse("1" %in% names(counts), counts["1"] * 1, 0)
  
  # Calculate allele frequencies
  total <- sum(counts)
  p0 <- allele0_count / (2 * total)
  p1 <- allele1_count / (2 * total)
  
  # Calculate the MAF
  maf <- min(p0, p1)
  
  # Check if MAF meets the threshold
  return(maf >= threshold)
}

# Apply MAF filter 1 to genotype matrix (row-wise)
maf_pass1 <- apply(geno_df, 1, sv_maf_filter1, threshold = 0.05)
filtered_geno_df1 <- geno_df[maf_pass1, ]

cat("Genotype matrix filtered with filter 1: retained", nrow(filtered_geno_df1), "variants\n")

# Save the filtered matrix from filter 1
write.table(filtered_geno_df1, file = "filtered_genotype_filter1.tsv", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# MAF Filter Function - 2nd Filter
sv_maf_filter2 <- function(geno_row, max_na = 10, threshold = 0.05) {
  # Count NA values in the row
  na_count <- sum(is.na(geno_row))
  
  # Discard variants with more than max_na NAs
  if (na_count > max_na) return(FALSE)
  
  # Replace remaining NAs with 0
  geno_row[is.na(geno_row)] <- 0
  
  # Count all alleles (0, 1, 2) — approximate allele counts
  counts <- table(geno_row)
  
  # Ensure counts for each genotype (0, 1, 2) are defined
  allele0_count <- ifelse("0" %in% names(counts), counts["0"] * 2, 0) +
    ifelse("1" %in% names(counts), counts["1"] * 1, 0)
  allele1_count <- ifelse("2" %in% names(counts), counts["2"] * 2, 0) +
    ifelse("1" %in% names(counts), counts["1"] * 1, 0)
  
  # Calculate allele frequencies
  total <- sum(counts)
  p0 <- allele0_count / (2 * total)
  p1 <- allele1_count / (2 * total)
  
  # Calculate the MAF
  maf <- min(p0, p1)
  
  # Check if MAF meets the threshold
  return(maf >= threshold)
}

# Apply MAF filter 2 to genotype matrix (row-wise)
maf_pass2 <- apply(geno_df, 1, sv_maf_filter2, threshold = 0.05)
filtered_geno_df2 <- geno_df[maf_pass2, ]

cat("Genotype matrix filtered with filter 2: retained", nrow(filtered_geno_df2), "variants\n")

# Save the filtered matrix from filter 2
write.table(filtered_geno_df2, file = "filtered_genotype_filter2.tsv", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Summary of results
cat("Filter 1 retained", nrow(filtered_geno_df1), "variants\n")
cat("Filter 2 retained", nrow(filtered_geno_df2), "variants\n")


# Expression filter: remove genes with low average expression
# You can adjust this threshold depending on your normalization method
expr_pass <- rowMeans(expr_df) > 1  # keep genes with mean VST > 1
filtered_expr_df <- expr_df[expr_pass, ]

cat("Expression matrix filtered: retained", nrow(filtered_expr_df), "genes\n")


# === Step 2: Convert filtered data.frames to SlicedData ===

# Write to temporary files so MatrixEQTL can load them
write.table(filtered_geno_df, "filtered_genotype.tsv", sep = "\t", quote = FALSE)
write.table(filtered_expr_df, "filtered_expression.tsv", sep = "\t", quote = FALSE)

# Function to load matrix into SlicedData
load_matrix <- function(file) {
  dat <- SlicedData$new()
  dat$fileDelimiter <- "\t"
  dat$fileOmitCharacters <- "NA"
  dat$fileSkipRows <- 1
  dat$fileSkipColumns <- 1
  dat$fileSliceSize <- 2000  # you can adjust this depending on memory
  dat$LoadFile(filename = file)
  return(dat)
}

# Load all 3 matrices
genotype_data <- load_matrix("cleaned_genotype_data.tsv")
expression_data <- load_matrix("cleaned_expression_data.tsv")
covariates_data <- load_matrix("cleaned_covariates_data.tsv")  # assumed to already be clean

# === Step 3: Run MatrixEQTL ===

useModel <- modelLINEAR  # You can also try modelANOVA or modelLINEAR_CROSS
output_file <- "MatrixEQTL_results.tsv"


#### cis-eQTLs only


# Load TIP and gene positions
TIP_positions <- read.table("svpos.tsv", header = TRUE, sep = "\t")
gene_positions <- read.table("genepos_with_isoform.tsv", header = TRUE, sep = "\t")

###load matrices
genotype_data <- load_matrix("cleaned_genotype_data.tsv")
expression_data <- load_matrix("cleaned_expression_data.tsv")
covariates_data <- load_matrix("cleaned_covariates_data.tsv") 


# Check format (must include SNP/gene IDs matching matrix row names, chr, and position)
# Example:
# snpspos.tsv:  snpid    chr    pos
# genepos.tsv:  geneid   chr    left    right


# Run only cis-eQTLs
cisDist <- 5e3  # 5kbp cis-window
cis_pv_threshold <- 1e-5

results <- Matrix_eQTL_main(
  snps = genotype_data,
  gene = expression_data,
  cvrt = covariates_data,
  output_file_name = output_file,
  pvOutputThreshold = 0,  # 0 = skip trans
  useModel = useModel,
  errorCovariance = numeric(),
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE,
  TIPspos = TIP_positions,
  genepos = gene_positions,
  cisDist = cisDist,
  pvOutputThreshold.cis = cis_pv_threshold
)

cis_eqtls <- results$cis$eqtls
write.table(cis_eqtls, file = "TIPs_cis_eQTL_results_new.tsv", sep = "\t", quote = FALSE, row.names = FALSE)



#### plotting any gene-TIP association


##### Load necessary libraries
library(ggplot2)

# File paths
expression_file <- "filtered_expression.tsv"
genotype_file <- "filtered_genotype.tsv"

# Load the filtered data correctly
expr_mat <- read.table(expression_file, header = TRUE, check.names = FALSE)
rownames(expr_mat) <- expr_mat[,1]  # Set the first column as row names
expr_mat <- expr_mat[,-1]  # Remove the first column after setting row names

geno_mat <- read.table(genotype_file, header = TRUE, row.names = 1, check.names = FALSE)

# Check the structure and first few rows after fixing
print("Expression matrix structure:")
print(str(expr_mat))
print(head(expr_mat))

print("Genotype matrix structure:")
print(str(geno_mat))
print(head(geno_mat))

# Function to plot expression (VST) per genotype
plot_eqtl <- function(gene_id, snp_id) {
  # Check if both gene and SNP exist in matrices
  if (gene_id %in% rownames(expr_mat) && snp_id %in% rownames(geno_mat)) {
    # Extract expression and genotype vectors
    expr_vec <- as.numeric(expr_mat[gene_id, ])
    geno_vec <- as.numeric(geno_mat[snp_id, ])
    
    # Print extracted vectors for debugging
    cat("Expression values for", gene_id, ":\n")
    print(expr_vec)
    cat("Genotype values for", snp_id, ":\n")
    print(geno_vec)
    
    # Combine into a data frame
    df <- data.frame(Genotype = factor(geno_vec), Expression = expr_vec)
    
    # Generate the boxplot
    p <- ggplot(df, aes(x = Genotype, y = Expression)) +
      geom_jitter(width = 0.2, size = 2, color = "steelblue") +
      geom_boxplot(alpha = 0.3, outlier.shape = NA) +
      labs(
        title = paste("Expression (VST) vs Genotype:", gene_id, "and", snp_id),
        x = "Genotype",
        y = "Expression (VST)"
      ) +
      theme_minimal()
    
    # Print and save the plot
    print(p)
    ggsave(paste0("eQTL_boxplot_", gene_id, "_vs_", snp_id, ".pdf"), plot = p, width = 6, height = 4)
    
  } else {
    cat("Either gene or SNP not found in matrices: ", gene_id, snp_id, "\n")
  }
}

# Example usage
plot_eqtl("Gh_A03G236800.1", "A03_109234090")

