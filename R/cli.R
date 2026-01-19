# ------------------------------------------------------------------
# cli.R
#
# Command-line interface for calculating duplex sequencing metrics.
#
# Defines the main CLI entrypoint, parses and validates command-line
# arguments, resolves input files, and delegates computation to
# `process_data()` in calculate.R.
#
# This script is executed via `main.R`.
# ------------------------------------------------------------------


source("R/calculate.R")


# Main CLI function
#' @export
main <- function() {
  # Parse arguments
  args <- parse_arguments()
  
  # Run metric computation
  res <- process_data(
    input     = args$input,
    output    = args$output,
    sample    = args$sample,
    rlen      = args$rlen,
    skips     = args$skips, 
    ref_fasta = args$ref_fasta, 
    metrics   = args$metrics, 
    cores     = args$cores
  )
  

  # Exit
  if (is.list(res) && isTRUE(res$success)) {
    cat("Process completed successfully\n")
    quit(status = 0)
  } else {
    cat("Process failed:", if (is.list(res)) res$error else "unknown error", "\n")
    quit(status = 1)
  }
}

# Parse command line arguments
parse_arguments <- function() {
  p <- argparse::ArgumentParser(description = "Calculate duplex QC metrics")
  
  p$add_argument("-i", "--input", nargs = "+", default = NULL,
                 help = "One or more input rinfo files (.txt or .txt.gz)")
  p$add_argument("--input_dir", default = "",
                 help = "Directory containing rinfo files")
  p$add_argument("--pattern", default = "\\.txt(\\.gz)?$",
                 help = "Regex pattern used with --input_dir (default: \\.txt(\\.gz)?$)")
  p$add_argument("--cores", type = "integer", default = 1,
                 help = "Number of cores for parallel processing (default: 1)")
  
  p$add_argument("-o", "--output", required = TRUE,
                 help = "Output CSV path (long format: sample,metric,value)")
  p$add_argument("-s", "--sample", default = NULL,
                 help = "Optional sample name(s). For multiple inputs, provide comma-separated names matching the number of files.")
  p$add_argument("--rlen", type = "integer", default = 151,
                 help = "Read length (default: 151)")
  p$add_argument("--skips", type = "integer", default = 5,
                 help = "Trimmed/ignored bases per read (NanoSeq = 5, xGen = 8)")
  p$add_argument("--ref_fasta", default = "",
                 help = "Optional reference genome FASTA. GC metrics are skipped if not provided.")
  p$add_argument(
    "--metrics", default = "",
    help = paste(
      "Comma-separated metric groups or metric names. Default: all.", 
      "Groups: gc,family.",
      "Individual: frac_singletons,efficiency,drop_out_rate.",
      "You may also specify grouped metrics (e.g. gc_single or family_mean) and the full group will be computed."
    )
  )

  
  args <- p$parse_args()
  
  
  if (!is.null(args$metrics) && nzchar(args$metrics)) {
    args$metrics <- gsub("\\s+", "", args$metrics)
  }
  
  # validate cores
  if (is.na(args$cores) || args$cores < 1) {
    stop("--cores must be >= 1")
  }
  
  # input handling 
  has_input     <- !is.null(args$input)
  has_input_dir <- nzchar(args$input_dir)
  
  if (has_input && has_input_dir) {
    stop("Specify either --input or --input_dir, not both")
  }
  if (!has_input && !has_input_dir) {
    stop("You must specify either --input or --input_dir")
  }
  
  if (has_input_dir) {
    if (!dir.exists(args$input_dir)) {
      stop("input_dir does not exist: ", args$input_dir)
    }
    files <- list.files(
      args$input_dir,
      pattern = args$pattern,
      full.names = TRUE
    )
    if (length(files) == 0) {
      stop("No input files found in input_dir matching pattern")
    }
    args$input <- files
  }
  
  # validate that all input files exist
  missing <- args$input[!file.exists(args$input)]
  if (length(missing) > 0) {
    stop("Input file(s) not found:\n", paste(missing, collapse = "\n"))
  }
  
  # normalise paths
  args$input <- normalizePath(args$input, mustWork = TRUE)

  
  # validate numeric CLI arguments early to avoid invalid metric computation
  if (is.na(args$rlen) || args$rlen <= 0)  stop("--rlen must be positive")
  if (is.na(args$skips) || args$skips < 0) stop("--skips must be >= 0")
  if (args$skips >= args$rlen) stop("--skips must be < --rlen")
  

  # output dir 
  odir <- dirname(args$output)
  if (!dir.exists(odir)) dir.create(odir, recursive = TRUE, showWarnings = FALSE)
  
  args
}

