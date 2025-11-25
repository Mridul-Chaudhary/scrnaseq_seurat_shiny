# install_packages.R
# Script to install all packages required for the Seurat-based single cell workflow.

# CRAN packages
cran_packages <- c(
  "Seurat",        # Core single-cell analysis
  "ggplot2",       # Visualization
  "devtools",      # For installing GitHub packages
  "usethis",       # Environment management
  "patchwork",     # Plot combination
  "remotes",       # Alternative for installing GitHub packages
  "scCATCH",       # Cell type annotation
  "openai",        # OpenAI API access
  "kable",         # Knit markdown tables
  "kableExtra",    # Knit markdown tables
  "dplyr"          # data manipulation
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c(
  "SingleR",       # Automated cell type annotation
  "celldex",       # Reference datasets for SingleR
  "scRNAseq"       # Reference datasets for scRNA-seq
)
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# GitHub packages
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Install 'presto' (fast marker detection)
if (!requireNamespace("presto", quietly = TRUE))
  devtools::install_github('immunogenomics/presto')

# Install 'Azimuth' (reference mapping)
if (!requireNamespace("Azimuth", quietly = TRUE))
  devtools::install_github("satijalab/azimuth", "seurat5")

# Install SeuratData from the correct branch
if (!requireNamespace("SeuratData", quietly = TRUE))
  devtools::install_github("satijalab/seurat-data", "seurat5")

# Install GPTCelltype (OpenAI-based annotation)
if (!requireNamespace("GPTCelltype", quietly = TRUE))
  remotes::install_github("Winnie09/GPTCelltype")

