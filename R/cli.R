# R/cli.R
suppressPackageStartupMessages({
  library(argparse)   # install.packages("argparse") or use renv::restore()
})

source("R/calculate.R")

# Main CLI function
#' @export
main <- function() {
  # Parse arguments
  args <- parse_arguments()

  # Setup logging
  # setup_logging(args$verbose)

  # Validate inputs
  # validate_inputs(args)
  
  # Call the worker
  res <- process_data(
    input     = args$input,
    output    = args$output,
    sample    = args$sample,
    rlen      = args$rlen,
    skips     = args$skips,
    ref_fasta = args$ref_fasta,
    skip_gc   = args$skip_gc,
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
                 help = "Sample ID (only valid when a single input file is used)")
  p$add_argument("--rlen", type = "integer", default = 151,
                 help = "Read length (default: 151)")
  p$add_argument("--skips", type = "integer", default = 5,
                 help = "Trimmed/ignored bases per read (NanoSeq = 5, xGen = 8)")
  p$add_argument("--ref_fasta", default = "",
                 help = "Optional reference genome FASTA (enables GC metrics)")
  p$add_argument("--skip_gc", default = "TRUE",
                 help = "TRUE/FALSE; disable GC even if FASTA is provided (default: TRUE)")
  p$add_argument("--metrics", default = "",
                 help = "Comma-separated metric groups or metric names (basic,gc,family). Default: all")
  p$add_argument("-v", "--verbose", action = "store_true",
                 help = "Verbose output")
  
  args <- p$parse_args()
  
  # normalise logicals / strings 
  args$skip_gc <- tolower(as.character(args$skip_gc)) %in%
    c("true", "1", "t", "yes", "y")
  
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
  
  # normalise paths (optional but recommended)
  args$input <- normalizePath(args$input, mustWork = TRUE)
  
  # sample sanity 
  if (length(args$input) > 1 && !is.null(args$sample) && nzchar(args$sample)) {
    stop("--sample can only be used with a single input file")
  }
  
  
  # output dir 
  odir <- dirname(args$output)
  if (!dir.exists(odir)) dir.create(odir, recursive = TRUE, showWarnings = FALSE)
  
  args
}

