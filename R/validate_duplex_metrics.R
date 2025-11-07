#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Compare CSV (efficiency, drop_out_rate) to ecoli_K12_metrics.rds reference
#
# Usage:
#    Rscript validate_duplex_metrics.R <csv> <ref_rds> <out_csv> [sample_id] [tolerance]
# ------------------------------------------------------------------

args <- commandArgs(TRUE)
if (length(args) < 3) stop("Usage: Rscript validate_duplex_metrics.R <csv> <ref_rds> <out_csv> [sample_id] [tolerance]")
csv_path <- args[1]; rds_path <- args[2]; out_path <- args[3]
sample_id <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NA_character_
tol <- if (length(args) >= 5 && nzchar(args[5])) as.numeric(args[5]) else 1e-3

norm <- function(x) {
  x <- tolower(gsub("[- ]+", "_", x))
  x <- sub("^value\\.", "", x)             # value.efficiency -> efficiency
  x <- ifelse(x %in% c("missed_fraction","dropout_rate"), "drop_out_rate", x)
  x <- ifelse(x %in% c("duplex_efficiency"), "efficiency", x)
  x
}

extract_pair <- function(df) {
  nms <- norm(names(df)); names(df) <- nms
  if (all(c("efficiency","drop_out_rate") %in% nms)) {
    return(c(efficiency=as.numeric(df[1,"efficiency"]), drop_out_rate=as.numeric(df[1,"drop_out_rate"])))
  }
  if (all(c("metric","value") %in% nms)) {
    mm <- norm(df$metric)
    return(c(
      efficiency   = as.numeric(df$value[match("efficiency", mm)]),
      drop_out_rate= as.numeric(df$value[match("drop_out_rate", mm, nomatch=match("dropout_rate", mm))])
    ))
  }
  stop("Unrecognised schema: ", paste(names(df), collapse=", "))
}

stopifnot(file.exists(csv_path), file.exists(rds_path))
our <- read.csv(csv_path, stringsAsFactors=FALSE)

if ("sample" %in% names(our)) {
  if (is.na(sample_id)) sample_id <- unique(our$sample)[1]
  our <- our[our$sample == sample_id, , drop=FALSE]
  if (nrow(our) == 0) stop("Sample not found in CSV: ", sample_id)
}
our_vals <- extract_pair(our)

ref <- readRDS(rds_path)
ref_df <- NULL
if (is.data.frame(ref)) ref_df <- ref else if (is.list(ref)) ref_df <- Filter(is.data.frame, ref)[[1]]
if (is.null(ref_df)) stop("Could not find a data.frame inside the RDS")

if ("sample" %in% names(ref_df) && !is.na(sample_id)) {
  ref_df <- ref_df[ref_df$sample == sample_id, , drop=FALSE]
  if (nrow(ref_df) == 0) stop("Sample not found in RDS: ", sample_id)
}
ref_vals <- extract_pair(ref_df)

cmp <- data.frame(
  sample   = sample_id,
  metric   = c("efficiency","drop_out_rate"),
  csv_val  = as.numeric(our_vals[c("efficiency","drop_out_rate")]),
  ref_val  = as.numeric(ref_vals[c("efficiency","drop_out_rate")]),
  abs_diff = NA_real_,
  pass     = NA,
  check.names = FALSE
)
cmp$abs_diff <- abs(cmp$csv_val - cmp$ref_val)
cmp$pass     <- cmp$abs_diff <= tol

dir.create(dirname(out_path), recursive=TRUE, showWarnings=FALSE)
write.csv(cmp, out_path, row.names=FALSE, quote=FALSE)

if (!all(cmp$pass, na.rm=TRUE)) {
  message("VALIDATION FAIL"); print(cmp); quit(status=1)
} else {
  message("VALIDATION PASS"); print(cmp); quit(status=0)
}

