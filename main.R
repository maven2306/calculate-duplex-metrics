#!/usr/bin/env Rscript

# Source required functions
source("R/cli.R")

# Install required packages if not available
required_packages <- c("argparse")
missing_packages <-
  required_packages[!required_packages %in% installed.packages()[, "Package"]]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:",
      paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cran.rstudio.com/")
}

# Load required libraries
suppressPackageStartupMessages({
  library(argparse)
})

# Run the main CLI function
main()
