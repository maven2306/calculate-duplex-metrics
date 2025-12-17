# ------------------------------------------------------------------
# calc_duplex_metrics.R
#
# Compute duplex-sequencing metrics from a single rinfo file and write
# them in long-format (sample, metric, value) CSV.
#
# This file is not a CLI entrypoint.
# It provides helper functions that are called from calculate.R
# (which is itself invoked via main.R -> cli.R).
#
# Description:
#   This module loads metric functions from efficiency_nanoseq_functions.R
#   and applies them to a single input rinfo file. All available metrics
#   returned by calculate_metrics_single() are exported, including:
#
#     - frac_singletons
#     - efficiency
#     - drop_out_rate
#     - GC metrics (gc_single, gc_both, gc_deviation)*
#     - family statistics (total_families, family_mean, family_median,
#       family_max, families_gt1, single_families, paired_families,
#       paired_and_gt1)
#
#   (*) GC metrics are only computed when a reference genome is provided.
#       If no reference is given (or GC is skipped), GC metrics are set to NA.
#
# Notes:
#   - Output is always written as long-format CSV suitable for MultiQC.
#   - Optional metric selection is handled by the calling code (calculate.R),
#     via the metrics_arg parameter.
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

# helper function that returns long-format data.frame (no writing)
calc_duplex_metrics_one_file_df <- function(
    input,
    sample = NULL,
    rlen = 151,
    skips = 5,
    ref_fasta = "",
    skip_gc = TRUE,
    metrics_arg = NULL,  # optional comma-separated metric names
    func_file = file.path("R", "efficiency_nanoseq_functions.R")
) {
  in_file <- normalizePath(input, mustWork = TRUE)
  
  if (is.null(sample) || !nzchar(sample)) {
    sample <- sub("\\.txt(\\.gz)?$", "", basename(in_file), TRUE)
  }
  
  rlen  <- as.integer(rlen)
  skips <- as.integer(skips)
  
  if (nzchar(ref_fasta)) {
    ref_fasta <- normalizePath(ref_fasta, mustWork = TRUE)
  }
  
  message("Using rlen=", rlen, " skips=", skips,
          " ref_fasta=", ifelse(nzchar(ref_fasta), ref_fasta, "<none>"))
  
  # removed rfunc_in: always source known file in repo
  if (!file.exists(func_file)) stop("Can't find efficiency_nanoseq_functions.R at: ", func_file)
  message("func_file: ", normalizePath(func_file))
  
  # load helpers into a private env that already has globals
  fn_env <- new.env(parent = .GlobalEnv)
  assign("rlen",  rlen,  envir = fn_env)
  assign("skips", skips, envir = fn_env)
  assign("bases_skipped", skips, envir = fn_env)
  sys.source(func_file, envir = fn_env)
  
  # if a reference FASTA was provided, expose it to the metrics code
  if (nzchar(ref_fasta)) {
    gf <- Rsamtools::FaFile(ref_fasta)
    
    fa <- Biostrings::readDNAStringSet(ref_fasta)
    fa_names <- sub("\\s.*$", "", names(fa))
    genome_max <- setNames(as.integer(Biostrings::width(fa)), fa_names)
    
    assign("genomeFile", gf,         envir = fn_env)
    assign("genome_max", genome_max, envir = fn_env)
  }
  
  # GC handling logic
  if (isTRUE(skip_gc) && exists("calculate_gc", envir = fn_env, inherits = FALSE)) {
    warning("GC calculation skipped by user flag — GC metrics will be set to NA.")
    fn_env$calculate_gc <- function(rbs, ...) {
      c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_)
    }
  } else if (!isTRUE(skip_gc) && !nzchar(ref_fasta) &&
             exists("calculate_gc", envir = fn_env, inherits = FALSE)) {
    stop("GC calculation requested but no reference genome provided. ",
         "Provide ref_fasta or set skip_gc=TRUE.")
  }
  
  # read rinfo
  rbs <- tryCatch(
    data.table::fread(in_file),
    error = function(e) e
  )
  if (inherits(rbs, "error")) stop("Failed to read input: ", rbs$message)
  
  # compute metrics (wide format, 1-row)
  tbl <- tryCatch(
    fn_env$calculate_metrics_single(rbs),
    error = function(e) e
  )
  if (inherits(tbl, "error")) stop("Metric calculation failed: ", tbl$message)
  
  if (is.null(tbl) || nrow(tbl) == 0) stop("Empty metrics table returned")
  
  metric_names <- names(tbl)
  
  # Optional metric selection (metrics_arg)
  if (!is.null(metrics_arg) && nzchar(metrics_arg)) {
    requested <- trimws(unlist(strsplit(metrics_arg, ",")))
    requested <- requested[nzchar(requested)]
    
    metric_names <- metric_names[metric_names %in% requested]
    
    if (length(metric_names) == 0) {
      stop("No valid metrics selected via metrics_arg")
    }
  }
  
  # return long-format df
  data.frame(
    sample = sample,
    metric = metric_names,
    value  = as.numeric(tbl[1, metric_names, drop = TRUE]),
    check.names = FALSE
  )
}

# call the df helper and write output 
calc_duplex_metrics_one_file <- function(
    input,
    output,
    sample = NULL,
    rlen = 151,
    skips = 5,
    ref_fasta = "",
    skip_gc = TRUE,
    metrics_arg = NULL,  # optional comma-separated metric names
    func_file = file.path("R", "efficiency_nanoseq_functions.R")
) {
  out_csv <- output
  
  out <- calc_duplex_metrics_one_file_df(
    input = input,
    sample = sample,
    rlen = rlen,
    skips = skips,
    ref_fasta = ref_fasta,
    skip_gc = skip_gc,
    metrics_arg = metrics_arg,
    func_file = func_file
  )
  
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  write.csv(out, out_csv, row.names = FALSE, quote = FALSE)
  cat("Wrote CSV:", out_csv, "\n")
  
  invisible(out_csv)
}

# compute many files, return one combined df 
calc_duplex_metrics_many_files_df <- function(
    inputs,
    rlen = 151,
    skips = 5,
    ref_fasta = "",
    skip_gc = TRUE,
    metrics_arg = NULL,
    cores = 1,
    func_file = file.path("R", "efficiency_nanoseq_functions.R")
) {
  if (is.null(inputs) || length(inputs) == 0) {
    stop("No input files provided")
  }
  
  inputs <- normalizePath(inputs, mustWork = TRUE)
  
  worker <- function(f) {
    calc_duplex_metrics_one_file_df(
      input = f,
      sample = NULL,          # default derived from filename
      rlen = rlen,
      skips = skips,
      ref_fasta = ref_fasta,
      skip_gc = skip_gc,
      metrics_arg = metrics_arg,
      func_file = func_file
    )
  }
  
  # macOS/Linux: use mclapply when cores > 1, Windows falls back to serial
  if (cores > 1 && .Platform$OS.type != "windows") {
    out_list <- parallel::mclapply(inputs, worker, mc.cores = cores)
  } else {
    out_list <- lapply(inputs, worker)
  }
  
  data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE) |> as.data.frame()
}

