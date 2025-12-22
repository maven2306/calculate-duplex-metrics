# ------------------------------------------------------------------
# calc_duplex_metrics.R
#
# Compute duplex-sequencing metrics from one or more rinfo files and
# return them in long-format (sample, metric, value)
#
# This file is not a CLI entrypoint.
# It provides helper functions that are called from calculate.R
# (which is invoked via main.R -> cli.R).
#
# Description:
#   This module loads metric functions from efficiency_nanoseq_functions.R
#   and applies them to rinfo files using a metric-group–based approach.
#
#   Metrics are organised into logical groups:
#     - individual   : frac_singletons, efficiency, drop_out_rate
#     - gc      : gc_single, gc_both, gc_deviation
#     - family  : total_families, family_mean, family_median, family_max,
#                 families_gt1, single_families, paired_families,
#                 paired_and_gt1
#
#   Metric selection is resolved before computation:
#     - If metrics_arg is NULL or empty, all metric groups are computed.
#     - If a group name is provided (e.g. "gc"), all metrics in that group
#       are computed together.
#     - If an individual metric name is provided (e.g. "gc_single"), it is
#       mapped to its corresponding group and the full group is computed.
#
#   Metrics are computed via calculate_metrics_selected(), avoiding the
#   "compute-all-then-filter" pattern.
#
# GC handling:
#   - GC metrics are computed only when the "gc" group is requested.
#   - A reference genome (ref_fasta) is required unless skip_gc = TRUE.
#   - If GC is skipped, GC metrics are returned as NA.
#
# Output:
#   - Results are returned in long format with columns:
#       sample, metric, value
# 
#
# Notes:
#   - This module performs no argument parsing.
#   - File I/O (CSV writing) and CLI interaction are handled upstream.
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


# call the df helper and write output 
calc_duplex_metrics_one_file_df <- function(
    input,
    sample = NULL,
    rlen = 151,
    skips = 5,
    ref_fasta = "",
    skip_gc = TRUE,
    metrics_arg = NULL,
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
  
  if (!file.exists(func_file)) stop("Can't find efficiency_nanoseq_functions.R at: ", func_file)
  message("func_file: ", normalizePath(func_file))
  
  fn_env <- new.env(parent = .GlobalEnv)
  assign("rlen",  rlen,  envir = fn_env)
  assign("skips", skips, envir = fn_env)
  assign("bases_skipped", skips, envir = fn_env)
  sys.source(func_file, envir = fn_env)
  

  sel <- fn_env$resolve_metric_selection(metrics_arg)
  groups_to_compute <- sel$groups
  individual_to_compute <- sel$individual
  
  
  # only load reference objects if GC metrics are requested
  if ("gc" %in% groups_to_compute && nzchar(ref_fasta)) {
    gf <- Rsamtools::FaFile(ref_fasta)
    
    fa <- Biostrings::readDNAStringSet(ref_fasta)
    fa_names <- sub("\\s.*$", "", names(fa))
    genome_max <- setNames(as.integer(Biostrings::width(fa)), fa_names)
    
    assign("genomeFile", gf,         envir = fn_env)
    assign("genome_max", genome_max, envir = fn_env)
  }
  
  # GC handling logic (only if GC is requested)
  if ("gc" %in% groups_to_compute && exists("calculate_gc", envir = fn_env, inherits = FALSE)) {
    if (isTRUE(skip_gc)) {
      warning("GC calculation skipped by user flag — GC metrics will be set to NA.")
      fn_env$calculate_gc <- function(rbs, ...) {
        c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_)
      }
    } else if (!nzchar(ref_fasta)) {
      stop("GC calculation requested but no reference genome provided. ",
           "Provide ref_fasta or set skip_gc=TRUE.")
    }
  }
  
  rbs <- tryCatch(
    data.table::fread(in_file),
    error = function(e) e
  )
  if (inherits(rbs, "error")) stop("Failed to read input: ", rbs$message)
  
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
  
  process_one_file <- function(f) {
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
    out_list <- parallel::mclapply(inputs, process_one_file, mc.cores = cores)
  } else {
    out_list <- lapply(inputs, process_one_file)
  }
  
  data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE) |> as.data.frame()
}

