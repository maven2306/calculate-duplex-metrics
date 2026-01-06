# ------------------------------------------------------------------
# calculate.R
#
# This script is the main orchestration layer for duplex sequencing
# metric calculation.
#
# It is invoked upstream via:
#   main.R -> cli.R -> calclate.R
#  
# It loads the core metric calculation functions, resolves which
# metric groups and individual metrics should be computed based on
# user input, and coordinates the processing of one or more input
# duplex data files.
#
# Key responsibilities:
#   - Loads metric calculation functions into a dedicated environment
#   - Resolves requested metric groups and individual metrics
#   - Handles optional GC metric computation using a reference genome
#   - Supports single-sample and multi-sample inputs
#   - Runs metric calculation serially or in parallel
#   - Writes results to a tidy CSV with one row per sample–metric pair
#
# All helper functions in this file are pure computation helpers:
# they do not perform argument parsing or command-line I/O. This
# script assumes inputs have already been validated upstream and
# focuses solely on coordinating metric computation and output.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(magrittr)    # %>%
  library(parallel)    # mclapply used upstream
  library(data.table)  # fread()
  library(R.utils)     # fread on .gz
  library(Biostrings)
  library(GenomicRanges)
  library(IRanges)
  library(Rsamtools)
  library(seqinr)
  suppressWarnings(try(library(dplyr), silent = TRUE)) # optional
})
`%>%` <- magrittr::`%>%`

# ------------------------------------------------------------------
# CHANGED - Pure helper functions: no argument parsing or file I/O
# ------------------------------------------------------------------

calc_duplex_metrics_one_file_df <- function(
    input,
    sample = NULL,
    rlen = 151,
    skips = 5,
    groups_to_compute,
    individual_to_compute,
    skip_gc = TRUE,
    ref_fasta = "",
    fn_env
) {
  # normalize path
  in_file <- normalizePath(input, mustWork = TRUE)
  
  if (is.null(sample) || !nzchar(sample)) {
    sample <- sub("\\.txt(\\.gz)?$", "", basename(in_file), TRUE)
  }
  
  # prepare reference genome only if GC metrics are requested
  if ("gc" %in% groups_to_compute && nzchar(ref_fasta)) {
    gf <- Rsamtools::FaFile(ref_fasta)
    fa <- Biostrings::readDNAStringSet(ref_fasta)
    fa_names <- sub("\\s.*$", "", names(fa))
    genome_max <- setNames(as.integer(Biostrings::width(fa)), fa_names)
    
    assign("genomeFile", gf, envir = fn_env)
    assign("genome_max", genome_max, envir = fn_env)
  }
  
  # CHANGED - GC override: if skip_gc, replace GC function with NA-returning stub
  if ("gc" %in% groups_to_compute && skip_gc && exists("calculate_gc", envir = fn_env, inherits = FALSE)) {
    fn_env$calculate_gc <- function(rbs, ...) {
      c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_)
    }
  }
  
  # read input
  rbs <- tryCatch(
    data.table::fread(in_file),
    error = function(e) e
  )
  if (inherits(rbs, "error")) stop("Failed to read input: ", rbs$message)
  
  # compute metrics
  tbl <- tryCatch(
    fn_env$calculate_metrics_selected(rbs, groups = groups_to_compute, individual = individual_to_compute),
    error = function(e) e
  )
  if (inherits(tbl, "error")) stop("Metric calculation failed: ", tbl$message)
  if (is.null(tbl) || nrow(tbl) == 0) stop("Empty metrics table returned")
  
  metric_names <- names(tbl)
  data.frame(
    sample = sample,
    metric = metric_names,
    value  = as.numeric(tbl[1, metric_names, drop = TRUE]),
    check.names = FALSE
  )
}


calc_duplex_metrics_many_files_df <- function(
    inputs,
    rlen = 151,
    skips = 5,
    groups_to_compute,
    individual_to_compute,
    skip_gc = TRUE,
    ref_fasta = "",
    cores = 1,
    fn_env
) {
  if (is.null(inputs) || length(inputs) == 0) stop("No input files provided")
  
  inputs <- normalizePath(inputs, mustWork = TRUE)
  
  process_one_file <- function(f) {
    calc_duplex_metrics_one_file_df(
      input = f,
      sample = NULL,
      rlen = rlen,
      skips = skips,
      groups_to_compute = groups_to_compute,
      individual_to_compute = individual_to_compute,
      skip_gc = skip_gc,
      ref_fasta = ref_fasta,
      fn_env = fn_env
    )
  }
  
  # CHANGED - parallelization logic unchanged
  if (cores > 1 && .Platform$OS.type != "windows") {
    out_list <- parallel::mclapply(inputs, process_one_file, mc.cores = cores)
  } else {
    out_list <- lapply(inputs, process_one_file)
  }
  
  data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE) |> as.data.frame()
}


# ------------------------------------------------------------------
# CHANGED - process_data now resolves metric selection, GC logic,
# and loads function environment. Helpers are pure and only compute metrics.
# ------------------------------------------------------------------

process_data <- function(
    input, output,
    sample = NULL,
    rlen = 151,
    skips = 5,
    ref_fasta = "",
    skip_gc = TRUE,
    metrics = "",
    cores = 1,
    func_file = file.path("R", "efficiency_nanoseq_functions.R")
) {
  # validate input
  if (is.null(input) || length(input) == 0) {
    return(list(success = FALSE, error = "No input files provided"))
  }
  
  if (is.na(cores) || cores < 1) {
    return(list(success = FALSE, error = "--cores must be >= 1"))
  }
  
  odir <- dirname(output)
  if (!dir.exists(odir)) dir.create(odir, recursive = TRUE, showWarnings = FALSE)
  
  # CHANGED - Load calculation functions once
  if (!file.exists(func_file)) {
    return(list(success = FALSE, error = paste0("Can't find efficiency_nanoseq_functions.R at: ", func_file)))
  }
  fn_env <- new.env(parent = .GlobalEnv)
  assign("rlen",  rlen,  envir = fn_env)
  assign("skips", skips, envir = fn_env)
  assign("bases_skipped", skips, envir = fn_env)
  sys.source(func_file, envir = fn_env)
  
  # CHANGED - Resolve metric selection once
  sel <- fn_env$resolve_metric_selection(metrics)
  groups_to_compute <- sel$groups
  individual_to_compute <- sel$individual
  
  # ---- single input ----
  if (length(input) == 1) {
    if (is.null(sample) || !nzchar(sample)) {
      sample <- sub("\\.txt(\\.gz)?$", "", basename(input), ignore.case = TRUE)
    }
    
    if (!file.exists(input)) {
      return(list(success = FALSE, error = paste0("Input not found: ", input)))
    }
    
    ok <- tryCatch({
      out_df <- calc_duplex_metrics_one_file_df(
        input = input,
        sample = sample,
        rlen = rlen,
        skips = skips,
        groups_to_compute = groups_to_compute,
        individual_to_compute = individual_to_compute,
        skip_gc = skip_gc,
        ref_fasta = ref_fasta,
        fn_env = fn_env
      )
      write.csv(out_df, output, row.names = FALSE, quote = FALSE)
      TRUE
    }, error = function(e) e)
    
    if (inherits(ok, "error")) return(list(success = FALSE, error = ok$message))
    
    return(list(success = TRUE))
  }
  
  # ---- multiple input ----
  if (!is.null(sample) && nzchar(sample)) {
    return(list(success = FALSE, error = "--sample can only be used with a single input file"))
  }
  
  missing <- input[!file.exists(input)]
  if (length(missing) > 0) {
    return(list(success = FALSE, error = paste0("Input file(s) not found:\n", paste(missing, collapse = "\n"))))
  }
  
  out_df <- tryCatch({
    calc_duplex_metrics_many_files_df(
      inputs = input,
      rlen = rlen,
      skips = skips,
      groups_to_compute = groups_to_compute,
      individual_to_compute = individual_to_compute,
      skip_gc = skip_gc,
      ref_fasta = ref_fasta,
      cores = cores,
      fn_env = fn_env
    )
  }, error = function(e) e)
  
  if (inherits(out_df, "error")) return(list(success = FALSE, error = out_df$message))
  
  ok <- tryCatch({
    write.csv(out_df, output, row.names = FALSE, quote = FALSE)
    TRUE
  }, error = function(e) e)
  
  if (inherits(ok, "error")) return(list(success = FALSE, error = paste0("Write failed: ", ok$message)))
  
  list(success = TRUE)
}
