# GIBGM Causality

Gene Interaction-Based Graphical Model for Causal Inference in Genomics

## Description

This project provides a novel method to infer causal structures among genes. 
It characterizes genes into causal genes and effect genes using functional data 
analysis (FDA) techniques, specifically Fourier expansion for genotype encoding.

## Background

Understanding the causal relationships between genes is crucial for uncovering 
the genetic architecture of complex diseases. Traditional methods often struggle 
with high-dimensional SNP (Single Nucleotide Polymorphism) data. This project 
addresses this challenge by:

1. Using Fourier basis expansion to encode genotype data
2. Applying PCA to determine the optimal number of basis functions
3. Inferring causal structures between gene pairs

## Features

- **Fourier Expansion Encoding**: Transforms SNP data using functional data analysis
- **PCA-based Basis Selection**: Automatically determines the number of Fourier basis functions
- **Gene-level Analysis**: Aggregates SNP-level information to gene-level representations

## Requirements

### R Version
- R >= 3.5.1

### R Packages
The following packages are required:

```r
# Core packages
library(fda)          # Functional data analysis
library(MASS)         # Statistical functions
library(stats)        # Basic statistical functions
library(factoextra)   # PCA visualization

# Optional packages for data handling
library(snpStats)     # SNP data handling
library(survival)     # Survival analysis
library(ggplot2)      # Plotting
```

## Installation

1. Install R (>= 3.5.1) from [CRAN](https://cran.r-project.org/)

2. Install required packages:

```r
install.packages(c("fda", "MASS", "factoextra", "ggplot2"))
```

3. For Bioconductor packages:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("snpStats")
```

## Project Structure

```
GIBGM_causality/
├── newFexpansion.R    # Main R script with Fourier expansion functions
└── README.md          # This file
```

## Usage

### Main Function: `fourierExpansion`

The core function performs Fourier expansion on genotype data:

```r
# Load required libraries
library(fda)
library(MASS)

# Perform Fourier expansion
result <- fourierExpansion(
  gene_idx = "GENE_NAME",    # Gene identifier
  geno = genotype_data,       # Genotype matrix
  gene_list = gene_annot,     # Gene annotation data
  snp_map = snp_map_data,     # SNP mapping information
  rng = 0                     # Range parameter
)
```

### Parameters

- `gene_idx`: Character string specifying the gene symbol to analyze
- `geno`: Data frame or matrix containing genotype data (samples x SNPs)
- `gene_list`: Data frame with gene annotation information
- `snp_map`: Data frame mapping SNPs to genomic positions
- `rng`: Numeric value for gene region extension (default: 0)

### Output

Returns a matrix with Fourier-expanded genotype values (zeta coefficients).

## Data Format

### Genotype Data
- Format: Data frame or matrix
- Rows: Samples/Individuals
- Columns: SNPs
- Values: Genotype coding (typically 0, 1, 2 for minor allele count)

### Gene Annotation
- `Gene_Symbol`: Gene identifiers
- `Chromosome`: Chromosome number
- `Start`: Gene start position
- `End`: Gene end position

### SNP Map
- `chr`: Chromosome
- `snp_id`: SNP identifier
- `posi`: Genomic position

## Methodology

1. **SNP Selection**: Select SNPs belonging to the target gene
2. **Position Normalization**: Normalize SNP positions to [0, 1] interval
3. **PCA Analysis**: Perform PCA on SNP data to determine variance explained
4. **Basis Selection**: Select number of basis functions (odd number based on 80% variance)
5. **Fourier Expansion**: Apply Fourier basis expansion to encode genotype data
6. **Output**: Return expansion coefficients (zeta)

## Example Workflow

```r
# Load libraries
library(fda)
library(MASS)

# Set data directory
dir_path <- '/path/to/your/data'

# Load data
geno_info <- read.table(file.path(dir_path, 'simGeno-chr2.raw'), header=TRUE)
gene_list <- read.csv(file.path(dir_path, 'gene.list.csv'))

# Process multiple genes
geno_expn <- matrix(0, nrow=nrow(snpm), ncol=0)
for(gene in unique(gene_list$Gene_Symbol)) {
  pergene <- fourierExpansion(gene, newgeno, gene_list, newmap, rng=0)
  colnames(pergene) <- paste('gene', gene, seq(1:ncol(pergene)), sep='_')
  geno_expn <- cbind(geno_expn, pergene)
  print(paste('Gene', gene, 'processed!'))
}
```

## Session Information

This project was developed and tested with:

```
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS 10.14

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
```

### Key Package Versions
- snpStats: 1.30.0
- survival: 2.42-6
- factoextra: 1.0.5
- ggplot2: 3.1.0
- fda: 2.4.8
- Matrix: 1.2-14
- MASS: 7.3-51

## References

1. Ramsay, J. O., & Silverman, B. W. (2005). Functional Data Analysis. Springer.
2. Storey, J. D., et al. (2005). Significance analysis of time course microarray experiments.

## License

MIT License

## Contact

For questions or issues, please open an issue on the repository.

## Acknowledgments

This work was inspired by the FRGEpistasis package and aims to improve 
gene-gene interaction analysis through functional data analysis techniques.
