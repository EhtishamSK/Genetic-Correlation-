# Genetic Correlation Analysis in R

## Overview
This R script performs **genetic correlation analysis** using **GBS-derived SNP data (HapMap format)** and **phenotypic trait data**. It is designed for plant breeders and geneticists who want to quantify the relationships between traits based on **additive genetic variance** estimated from the **genomic relationship matrix (GRM)**.

The script allows users to:

- Import and process **GBS SNP data** in HapMap format.
- Convert genotype calls (`AA`, `AG`, `GG`) into **numeric values** (0, 1, 2).
- Construct the **genomic relationship matrix (GRM)** using the `rrBLUP` package.
- Fit **multivariate mixed models** using `sommer`.
- Estimate **variance components** and **genetic correlations** between traits.
- Generate a **genetic correlation heatmap** for visualization.
- Export results as `.csv` files for reporting and downstream use.

This workflow enables breeders to **understand trait-trait relationships** at the genetic level, helping in **selection strategies** and **multi-trait improvement**.

---

## Features

- Supports **multiple traits** in one dataset.
- Fits **multivariate models** with genomic relationship matrix (GRM).
- Estimates **variance-covariance components**.
- Calculates **genetic correlations** between all trait pairs.
- Saves results as:
  - **Pairwise correlations** (`genetic_correlations.csv`)
  - **Correlation matrix** (`genetic_correlation_matrix.csv`)
  - **Correlation matrix with significance** (`genetic_correlations_with_significance.csv`)
- Generates a **heatmap plot** of genetic correlations.

---

## Input Requirements

1. **Genotypic Data (HapMap format)**  
   - SNP file in HapMap format (`.hmp.txt`)  
   - The script skips the first 10 header lines automatically.  
   - Genotype calls: `AA`, `AG/GA`, `GG`.

2. **Phenotypic Data (CSV file)**  
   - Must include:
     - `geno` – genotype names (matching HapMap file)  
     - Trait columns (numeric values for measured traits)  

> ⚠️ **Important:** Genotype names in the phenotypic file must match the HapMap file to align correctly with the GRM.

---

## Output

The script generates:

- **Genetic correlation values** for all trait pairs (`genetic_correlations.csv`).
- **Genetic correlation matrix** across all traits (`genetic_correlation_matrix.csv`).
- **Heatmap plot** of genetic correlations (`ggplot2` visualization).
- **Optional significance testing** with p-values and annotated correlation matrix.


---

## Citation

If you use this script, please cite the underlying R packages:

- Endelman, J.B. (2011). **Ridge Regression and Other Kernels for Genomic Selection with R Package rrBLUP.** *Plant Genome* 4(3):250-255.  
[https://doi.org/10.3835/plantgenome2011.08.0024](https://doi.org/10.3835/plantgenome2011.08.0024)

- Covarrubias-Pazaran, G. (2016). **Genome-Assisted Prediction of Quantitative Traits Using the R Package sommer.** *PLoS ONE* 11(6): e0156744.  
[https://doi.org/10.1371/journal.pone.0156744](https://doi.org/10.1371/journal.pone.0156744)

---

## Author

**Ehtisham Khokhar**  
New Mexico State University  
Email: ehtishamshakeel@gmail.com

