#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# calc_duplex_metrics.R
#
# Compute duplex-sequencing metrics from a single rinfo file and write
# them in long-format (sample, metric, value) CSV.
#
# Usage:
#     Rscript calc_duplex_metrics.R \
#     <input_rinfo> <output.csv> \
#     [sample_id] [rfunc_dir] [rlen] [skips] [ref_fasta] [skip_gc]
#
# Description:
#   This script loads metric functions from efficiency_nanoseq_functions.R
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
#   - This script supports both wide (column-based) and long (metric/value)
#     schemas returned by upstream functions.
#   - Output is always written as long format CSV suitable for MultiQC.
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

# 1. Parse CLI args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(paste(
    "Usage:",
    "Rscript calc_duplex_metrics.R <input_rinfo.txt(.gz)> <output.csv> [sample_id] [rfunc_dir_or_file] [rlen] [skips] [ref_fasta] [skip_gc]",
    sep = "\n"
  ))
}

in_file  <- normalizePath(args[1], mustWork = TRUE)
out_csv  <- args[2]
sample   <- if (length(args) >= 3 && nzchar(args[3])) args[3] else sub("\\.txt(\\.gz)?$", "", basename(in_file), TRUE)
rfunc_in <- if (length(args) >= 4 && nzchar(args[4])) args[4] else Sys.getenv("DUPLEX_RFUNCDIR", unset = "")


rlen  <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 151
skips <- if (length(args) >= 6 && nzchar(args[6])) as.integer(args[6]) else 5

# Reference FASTA
ref_fasta <- if (length(args) >= 7 && nzchar(args[7])) normalizePath(args[7]) else ""

# New CLI flag: -- skip_gc
# If --skip_gc = TRUE -> don’t run GC even if reference is available
# If reference is missing AND user did not skip GC -> show error
skip_gc <- FALSE
if (length(args) >= 8) skip_gc <- args[8] %in% c("TRUE", "true", "1")


message("Using rlen=", rlen, " skips=", skips,
        " ref_fasta=", ifelse(nzchar(ref_fasta), ref_fasta, "<none>"))


# Locate function file
func_file <- if (nzchar(rfunc_in) && grepl("\\.R$", rfunc_in, ignore.case = TRUE)) rfunc_in else file.path(rfunc_in, "efficiency_nanoseq_functions.R")
if (!file.exists(func_file)) stop("Can't find efficiency_nanoseq_functions.R at: ", func_file)
message("func_file: ", normalizePath(func_file))


# Load helpers into a private env that already has globals
fn_env <- new.env(parent = .GlobalEnv)
assign("rlen",  rlen,  envir = fn_env)
assign("skips", skips, envir = fn_env)
assign("bases_skipped", skips, envir = fn_env)
sys.source(func_file, envir = fn_env)

# If a reference FASTA was provided, expose it to the metrics code
if (nzchar(ref_fasta)) {
  # open FASTA as a FaFile for scanFa()
  gf <- Rsamtools::FaFile(ref_fasta)
  
  # derive chromosome lengths (genome_max) from the FASTA
  fa <- Biostrings::readDNAStringSet(ref_fasta)
  # Take only the first whitespace-delimited token from each FASTA header
  fa_names <- sub("\\s.*$", "", names(fa))
  genome_max <- setNames(as.integer(Biostrings::width(fa)), fa_names)
  
  assign("genomeFile", gf,        envir = fn_env)
  assign("genome_max", genome_max, envir = fn_env)
}



# MVP
# GC handling logic:
# - if skip_gc = TRUE,force GC metrics to NA, regardless of reference
# - if skip_gc = FALSE and ref_fasta provided, run real GC (calculate_gc as defined in functions file)
# - if skip_gc = FALSE and NO ref_fasta, error (user asked for GC but gave no reference)
if (skip_gc && exists("calculate_gc", envir = fn_env, inherits = FALSE)) {
  warning("GC calculation skipped by user flag — GC metrics will be set to NA.")
  fn_env$calculate_gc <- function(rbs, ...) {
    c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_)
  }
} else if (!skip_gc && !nzchar(ref_fasta) && exists("calculate_gc", envir = fn_env, inherits = FALSE)) {
  stop("GC calculation requested but no reference genome provided. ",
       "Provide ref_fasta or set skip_gc=TRUE.")
}

# 4. Compute metrics for exactly this file (match on basename w/o .txt/.gz)
pattern <- sub("\\.txt(\\.gz)?$", "", basename(in_file), TRUE)
res <- fn_env$calc_metrics_new_rbs(rinfo_dir = dirname(in_file), pattern = pattern, cores = 1)
if (is.null(res)) stop("No result returned.")

# calc_metrics_new_rbs (as written) returns a list (one element per matched file)
tbl <- if (is.list(res) && !is.data.frame(res)) res[[1]] else res
if (is.null(tbl) || nrow(tbl) == 0) stop("Empty metrics table.")

# 5. Extract Efficiency & Drop_out_rate (supports wide or long schemas)
#if (all(c("efficiency","drop_out_rate") %in% names(tbl))) {
#eff <- as.numeric(tbl$efficiency[1]); dr <- as.numeric(tbl$drop_out_rate[1])
#} else if (all(c("metric","value") %in% names(tbl))) {
#  mm <- tolower(gsub("[- ]","_", tbl$metric))
#  eff <- as.numeric(tbl$value[match("efficiency", mm)])
#  dr  <- as.numeric(tbl$value[match("drop_out_rate", mm, nomatch = match("dropout_rate", mm))])
#} else {
#  stop("Unrecognised metrics table schema: ", paste(names(tbl), collapse = ", "))
#}


# Include more metrics from efficiency_nanoseq_functons.R
metrics <- c(
  "frac_singletons",
  "efficiency",
  "drop_out_rate",
  "gc_single",
  "gc_both",
  "gc_deviation",
  "total_families",
  "family_mean",
  "family_median",
  "family_max",
  "families_gt1",
  "single_families",
  "paired_families",
  "paired_and_gt1"
            )

# Extract metrics (supports wide or long schemas)
if (any(metrics %in% names(tbl))) {
  metric_values <- setNames(
    as.numeric(tbl[1, metrics[metrics %in% names(tbl)], drop = TRUE]),
    metrics[metrics %in% names(tbl)]
  )
} else if (all(c("metric", "value") %in% names(tbl))) {
  mm <- tolower(gsub("[- ]","_", tbl$metric))
  metric_values <- vapply(
    metrics,
    function(m) {
      idx <- match(m, mm)
      if (is.na(idx)) NA_real_ else as.numeric(tbl$value[idx])
    },
    numeric(1)
  )
  names(metric_values) <- metrics
  
} else {
  stop("Unrecognised metrics table schema: ", paste(names(tbl), collapse = ", "))
}


# Write long-format CSV (sample, metric, value)
out <- data.frame(sample = sample,
                  #metric = c("efficiency","drop_out_rate"),
                  metric = names(metric_values),
                  #value  = c(eff, dr),
                  value = as.numeric(metric_values),
                  check.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(out, out_csv, row.names = FALSE, quote = FALSE)
cat("Wrote CSV:", out_csv, "\n")
