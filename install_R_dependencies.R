#!/usr/bin/env Rscript
# Install R dependencies required to reproduce the analyses.
# Tested on R >= 4.3.0 (CRAN snapshot 2026-04-15).

cran_packages <- c(
  # Data wrangling
  "dplyr", "tidyr", "stringr", "purrr", "readr",
  # Survival / causal inference
  "survival", "survey", "rms", "MatchIt", "WeightIt", "cobalt",
  "PSweight", "cmprsk", "boot",
  # Modeling / regularization
  "glmnet",
  # Visualization
  "ggplot2", "gridExtra", "ggpubr", "patchwork", "scales",
  # Reporting / tables
  "gtsummary", "officer", "flextable",
  # I/O
  "openxlsx", "data.table"
)

installed <- installed.packages()[, "Package"]
to_install <- setdiff(cran_packages, installed)

if (length(to_install) > 0) {
  cat("Installing:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://cloud.r-project.org/")
} else {
  cat("All required R packages are already installed.\n")
}
