#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# calc_duplex_metrics.R
#
# Compute duplex-sequencing metrics from a single rinfo file and write
# them in long-format (sample, metric, value) CSV.
#
# Usage:
#   Rscript calc_duplex_metrics.R \
#     <input_rinfo.txt(.gz)> <output.csv> \
#     [sample_id] [rfunc_dir_or_file] [rlen] [skips]
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
#   (*) GC metrics are returned as NA at the MVP stage unless a reference
#       genome and required Bioconductor dependencies are provided.
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
  suppressWarnings(try(library(dplyr), silent = TRUE)) # optional
})
`%>%` <- magrittr::`%>%`

# 1. Parse CLI args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(paste(
    "Usage:",
    "Rscript calc_duplex_metrics.R <input_rinfo.txt(.gz)> <output.csv> [sample_id] [rfunc_dir_or_file] [rlen] [skips]",
    sep = "\n"
  ))
}
in_file  <- normalizePath(args[1], mustWork = TRUE)
out_csv  <- args[2]
sample   <- if (length(args) >= 3 && nzchar(args[3])) args[3] else sub("\\.txt(\\.gz)?$", "", basename(in_file), TRUE)
rfunc_in <- if (length(args) >= 4 && nzchar(args[4])) args[4] else Sys.getenv("DUPLEX_RFUNCDIR", unset = "")

rlen  <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 151
skips <- if (length(args) >= 6 && nzchar(args[6])) as.integer(args[6]) else 5
message("Using rlen=", rlen, " skips=", skips)

# 2. Locate function file
func_file <- if (nzchar(rfunc_in) && grepl("\\.R$", rfunc_in, ignore.case = TRUE)) rfunc_in else file.path(rfunc_in, "efficiency_nanoseq_functions.R")
if (!file.exists(func_file)) stop("Can't find efficiency_nanoseq_functions.R at: ", func_file)
message("func_file: ", func_file)

# 3. Load helpers into a private env that already has globals
fn_env <- new.env(parent = .GlobalEnv)
assign("rlen",  rlen,  envir = fn_env)
assign("skips", skips, envir = fn_env)
assign("bases_skipped", skips, envir = fn_env)
sys.source(func_file, envir = fn_env)

# MVP: GC not needed; neutralise to avoid genome_* / Bioconductor deps 
if (exists("calculate_gc", envir = fn_env, inherits = FALSE)) {
  fn_env$calculate_gc <- function(rbs, ...) {
    c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_)
  }
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

