#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# main.R
#
# CLI entrypoint for the duplex metric calculation tool.
#
# This script initialises the runtime environment, ensures required
# packages are available, sources the CLI implementation, and invokes
# the `main()` function defined in cli.R.
# ------------------------------------------------------------------

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
