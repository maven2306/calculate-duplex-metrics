# R/cli.R
suppressPackageStartupMessages({
  library(argparse)   # install.packages("argparse") or use renv::restore()
})

source("R/calculate.R")

#' Main CLI function
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
    rfunc_dir = args$rfunc_dir,
    rlen      = args$rlen,
    skips     = args$skips
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

#' Parse command line arguments
parse_arguments <- function() {
  p <- argparse::ArgumentParser(description = "Calculate duplex QC metrics")
  
  p$add_argument("-i", "--input",  required = TRUE,
                 help = "Input rinfo file (.txt or .txt.gz)")
  p$add_argument("-o", "--output", required = TRUE,
                 help = "Output CSV path (long format: sample,metric,value)")
  p$add_argument("-s", "--sample", default = NULL,
                 help = "Sample ID (defaults to input basename without .txt[.gz])")
  p$add_argument("--rfunc_dir", required = TRUE,
                 help = "Path to G000204_duplex/code OR full path to efficiency_nanoseq_functions.R")
  p$add_argument("--rlen",  type = "integer", default = 151,
                 help = "Read length (default: 151)")
  p$add_argument("--skips", type = "integer", default = 5,
                 help = "Trimmed/ignored bases per read (Nano=5, XGen=8)")
  p$add_argument("-v", "--verbose", action = "store_true",
                 help = "Verbose output")
  args <- p$parse_args()
  
# ---- derive defaults / validate ----
  if (is.null(args$sample) || !nzchar(args$sample)) {
    bn <- basename(args$input)
    args$sample <- sub("\\.txt(\\.gz)?$", "", bn, ignore.case = TRUE)
  }
  
  if (!file.exists(args$input)) {
    stop("Input file not found: ", args$input)
  }
  
  odir <- dirname(args$output)
  if (!dir.exists(odir)) dir.create(odir, recursive = TRUE, showWarnings = FALSE)
  
  # rfunc_dir can be a directory OR the exact .R file
  if (!file.exists(args$rfunc_dir)) {
    stop("--rfunc_dir path does not exist: ", args$rfunc_dir)
  }
  if (dir.exists(args$rfunc_dir)) {
    fn <- file.path(args$rfunc_dir, "efficiency_nanoseq_functions.R")
    if (!file.exists(fn)) {
      stop("efficiency_nanoseq_functions.R not found in --rfunc_dir: ", args$rfunc_dir)
    }
  }
  
  # basic numeric sanity
  if (is.na(args$rlen) || args$rlen <= 0)  stop("--rlen must be positive")
  if (is.na(args$skips) || args$skips < 0) stop("--skips must be >= 0")
  
  args
}
