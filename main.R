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


# Check required packages are available (install via renv::restore())
required_packages <- c("argparse")
missing_packages <-
  required_packages[!required_packages %in% installed.packages()[, "Package"]]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "), "\n",
    "Please run renv::restore() before running this script.",
    call. = FALSE
  )
}

# Load required libraries
suppressPackageStartupMessages({
  library(argparse)
})

# Source CLI and run 
source("R/cli.R")
main()
