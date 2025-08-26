#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(magrittr)    # %>%
  library(parallel)    # mclapply used upstream
  library(data.table)  # fread()
  library(R.utils)     # fread on .gz
  suppressWarnings(try(library(dplyr), silent = TRUE)) # optional
})
`%>%` <- magrittr::`%>%`

# ---- args ----
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

# ---- resolve functions file ----
func_file <- if (nzchar(rfunc_in) && grepl("\\.R$", rfunc_in, ignore.case = TRUE)) rfunc_in else file.path(rfunc_in, "efficiency_nanoseq_functions.R")
if (!file.exists(func_file)) stop("Can't find efficiency_nanoseq_functions.R at: ", func_file)
message("func_file: ", func_file)

# ---- source upstream into env that already has globals ----
fn_env <- new.env(parent = .GlobalEnv)
assign("rlen",  rlen,  envir = fn_env)
assign("skips", skips, envir = fn_env)
assign("bases_skipped", skips, envir = fn_env)
sys.source(func_file, envir = fn_env)

# MVP: GC not needed; neutralise to avoid genome_* / Bioconductor deps
if (exists("calculate_gc", envir = fn_env, inherits = FALSE)) {
  fn_env$calculate_gc <- function(rbs, ...) NA_real_
}

# ---- run on exactly this file ----
pattern <- sub("\\.txt(\\.gz)?$", "", basename(in_file), TRUE)
res <- fn_env$calc_metrics_new_rbs(rinfo_dir = dirname(in_file), pattern = pattern, cores = 1)
if (is.null(res)) stop("No result returned.")

# calc_metrics_new_rbs (as written) returns a list (one element per matched file)
tbl <- if (is.list(res) && !is.data.frame(res)) res[[1]] else res
if (is.null(tbl) || nrow(tbl) == 0) stop("Empty metrics table.")

# ---- extract efficiency & drop_out_rate (wide or long) ----
if (all(c("efficiency","drop_out_rate") %in% names(tbl))) {
  eff <- as.numeric(tbl$efficiency[1]); dr <- as.numeric(tbl$drop_out_rate[1])
} else if (all(c("metric","value") %in% names(tbl))) {
  mm <- tolower(gsub("[- ]","_", tbl$metric))
  eff <- as.numeric(tbl$value[match("efficiency", mm)])
  dr  <- as.numeric(tbl$value[match("drop_out_rate", mm, nomatch = match("dropout_rate", mm))])
} else {
  stop("Unrecognised metrics table schema: ", paste(names(tbl), collapse = ", "))
}

# ---- write long-format CSV ----
out <- data.frame(sample = sample,
                  metric = c("efficiency","drop_out_rate"),
                  value  = c(eff, dr),
                  check.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(out, out_csv, row.names = FALSE, quote = FALSE)
cat("Wrote CSV:", out_csv, "\n")

