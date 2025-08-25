# ------------------------------------------------------------
# Genetic Correlation Analysis using GBS Data and Phenotypic Data
# Author: Ehtisham Khokhar
# Email: ehtisham@example.com
#
# Description:
# This script performs genetic correlation analysis using SNP data 
# (HapMap format from GBS) and phenotypic trait data. It:
#   1. Processes HapMap SNP data
#   2. Converts genotypes into numeric format
#   3. Builds the Genomic Relationship Matrix (GRM)
#   4. Fits multivariate mixed models using sommer
#   5. Estimates genetic correlations among traits
#   6. Exports results and generates a heatmap
# ------------------------------------------------------------

# ------------------------
# 1. Set working directory
# ------------------------
setwd("C:/Users/ehtis/OneDrive - New Mexico State University/SUNNY/Research Projects/Mechanical Harvest Projects/Phenotype manuscript/single location/NSH")
getwd()

# ---------------------------------------------------------
# 2. Load and process HapMap SNP data (skip first 10 lines)
# ---------------------------------------------------------
# My GBS file is in HapMap version, therefore we skip first 10 lines
hapmap_data <- read.table("GBS_PANEL_DIPLOID_75G.hmp.txt", header = FALSE, skip = 10)

# Quick check of data
head(hapmap_data)
str(hapmap_data)

# Load library for data handling
library(dplyr)

# Extract genotype columns (from V12 onwards in HapMap)
genotype_data <- hapmap_data[, 12:ncol(hapmap_data)]
str(genotype_data)

# ---------------------------------------------------------
# 3. Convert genotype codes (AA=0, AG=1, GG=2)
# ---------------------------------------------------------
convert_genotype <- function(genotype) {
  if (is.na(genotype)) {
    return(NA)  
  } else if (genotype == "AA") {
    return(0)
  } else if (genotype == "AG" | genotype == "GA") {
    return(1)
  } else if (genotype == "GG") {
    return(2)
  } else {
    return(NA)  
  }
}

# Apply conversion function across genotype dataset
numeric_genotypes <- apply(genotype_data, MARGIN = c(1, 2), FUN = convert_genotype)

# Transpose: rows = genotypes, columns = SNPs
numeric_genotypes <- t(numeric_genotypes)

# ---------------------------------------------------------
# 4. Build Genomic Relationship Matrix (GRM)
# ---------------------------------------------------------
library(rrBLUP)

# A.mat expects rows=genotypes and cols=SNPs
G <- A.mat(numeric_genotypes, return.impute = TRUE)  
str(G)
str(G$A)  

# Extract genotype names from GRM
geno_names_grm <- rownames(G$A)

# ---------------------------------------------------------
# 5. Import Phenotypic Data
# ---------------------------------------------------------
pheno <- read.csv("mydata_gc.csv")
str(pheno)
pheno$geno <- as.factor(pheno$geno)

# Extract genotype names
geno_names_pheno <- pheno$geno  

# Check alignment with GRM
all(geno_names_pheno == geno_names_grm)

# If mismatched, reorder phenotypic data (uncomment if needed)
# pheno <- pheno[match(geno_names_grm, pheno$geno), ]

# ---------------------------------------------------------
# 6. Fit Multivariate Models using sommer
# ---------------------------------------------------------
library(sommer)

# Example model with two traits (PWDT, PHT)
model_multivariate1 <- mmer(cbind(PWDT, PHT) ~ 1,
                            random = ~ vsr(geno, Gu = G$A),
                            rcov = ~ units,
                            data = pheno)

# Alternative with tolerance parameter
model_multivariate2 <- mmer(cbind(PWDT, PHT) ~ 1,
                            random = ~ vsr(geno, Gu = G$A),
                            rcov = ~ units,
                            tolParInv = 1e-2,
                            data = pheno)

# Extract variance-covariance components
varcov_components <- summary(model_multivariate2)$varcomp

# Example: covariance and variances for traits
cov_A_PWDT_PHT <- varcov_components["u:geno.PWDT-PHT", "VarComp"]
sigma_A_PWDT <- varcov_components["u:geno.PWDT-PWDT", "VarComp"]
sigma_A_PHT <- varcov_components["u:geno.PHT-PHT", "VarComp"]

# Compute genetic correlation
genetic_correlation_PWDT_PHT <- cov_A_PWDT_PHT / sqrt(sigma_A_PWDT * sigma_A_PHT)
genetic_correlation_PWDT_PHT

# ---------------------------------------------------------
# 7. Automate Genetic Correlation for All Trait Pairs
# ---------------------------------------------------------
genetic_correlation_list <- list()
traits <- colnames(pheno)[2:8]   # adjust range based on your data

# Loop over all trait pairs
for (i in 1:(length(traits)-1)) {
  for (j in (i+1):length(traits)) {
    trait1 <- traits[i]
    trait2 <- traits[j]
    
    model_multivariate <- tryCatch({
      mmer(as.formula(paste("cbind(", trait1, ",", trait2, ") ~ 1")),
           random = ~ vsr(geno, Gu = G$A),
           rcov = ~ units,
           data = pheno)
    }, error = function(e) {
      message(paste("Error in fitting model for", trait1, "and", trait2, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(model_multivariate)) {
      varcov_components <- summary(model_multivariate)$varcomp
      if (paste0("u:geno.", trait1, "-", trait2) %in% rownames(varcov_components)) {
        cov_A <- varcov_components[paste0("u:geno.", trait1, "-", trait2), "VarComp"]
        sigma1 <- varcov_components[paste0("u:geno.", trait1, "-", trait1), "VarComp"]
        sigma2 <- varcov_components[paste0("u:geno.", trait2, "-", trait2), "VarComp"]
        genetic_correlation <- cov_A / sqrt(sigma1 * sigma2)
        genetic_correlation_list[[paste(trait1, trait2, sep = "_")]] <- genetic_correlation
      }
    }
  }
}

# Save results
genetic_correlation_df <- as.data.frame(do.call(rbind, genetic_correlation_list))
colnames(genetic_correlation_df) <- "Genetic_Correlation"
write.csv(genetic_correlation_df, "genetic_correlations.csv", row.names = TRUE)

# ---------------------------------------------------------
# 8. Build and Save Correlation Matrix
# ---------------------------------------------------------
# (Code for building correlation matrix, filling upper/lower triangle, saving as CSV)
# ---------------------------------------------------------

# ---------------------------------------------------------
# 9. Visualize Genetic Correlations (Heatmap)
# ---------------------------------------------------------
library(ggplot2)
library(reshape2)

# Keep upper triangle for plotting
upper_triangle <- function(correlation_matrix) {
  correlation_matrix[lower.tri(correlation_matrix)] <- NA
  return(correlation_matrix)
}

upper_corr_matrix <- upper_triangle(as.matrix(correlation_matrix_df))
genetic_correlation_long <- melt(upper_corr_matrix, na.rm = TRUE)

# Heatmap with values displayed
heatmap_plot <- ggplot(data = genetic_correlation_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, space = "Lab", 
                       name="Genetic Correlations") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  labs(x = "Traits", y = "Traits", title = "Genetic Correlation Heatmap") +
  coord_fixed()

print(heatmap_plot)

# ---------------------------------------------------------
# 10. Extension: Add Significance Testing for Genetic Correlations
# ---------------------------------------------------------
# (Code for calculating p-values, significance, and matrix with significance marks)
# ---------------------------------------------------------

