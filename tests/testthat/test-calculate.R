library(testthat)

# ----------------------
# Metrics to validate
# ---------------------
metrics <- c(
  "frac_singletons", "efficiency", "drop_out_rate","total_families", "family_mean",
  "family_median", "family_max", "families_gt1", "single_families",
  "paired_families", "paired_and_gt1"
)

# Normalize column/metric names
norm <- function(x) {
  x <- tolower(gsub("[- ]+", "_", x))
  x <- sub("^value\\.", "", x)
  x <- ifelse(x %in% c("missed_fraction", "dropout_rate"), "drop_out_rate", x)
  ifelse(x %in% "duplex_efficiency", "efficiency", x)
}

# Extract metrics from a data frame (wide or long)
extract_metrics <- function(df, metrics) {
  names(df) <- norm(names(df))

  if (all(metrics %in% names(df))) {
    return(setNames(as.numeric(df[1, metrics]), metrics))
  }

  if (all(c("metric", "value") %in% names(df))) {
    mm <- norm(df$metric)
    vals <- sapply(metrics, function(m) as.numeric(df$value[match(m, mm)]))
    return(vals)
  }

  stop("Unrecognised schema: ", paste(names(df), collapse = ", "))
}

# Validate samples
validate_metrics <- function(csv_path, rds_path, metrics, tol = 1e-3) {
  stopifnot(file.exists(csv_path), file.exists(rds_path))

  csv_df <- read.csv(csv_path, stringsAsFactors = FALSE)
  rds_raw <- readRDS(rds_path)
  rds_df <- if (is.data.frame(rds_raw)) rds_raw else Filter(is.data.frame, rds_raw)[[1]]

  if (!"sample" %in% names(csv_df)) stop("CSV missing 'sample' column")
  if (!"sample" %in% names(rds_df)) stop("RDS missing 'sample' column")

  samples <- intersect(unique(csv_df$sample), unique(rds_df$sample))
  if (length(samples) == 0) stop("No common samples between CSV and RDS")

  do.call(rbind, lapply(samples, function(s) {
    csv_sub <- csv_df[csv_df$sample == s, , drop = FALSE]
    rds_sub <- rds_df[rds_df$sample == s, , drop = FALSE]

    csv_vals <- extract_metrics(csv_sub, metrics)
    rds_vals <- extract_metrics(rds_sub, metrics)

    data.frame(
      sample = s,
      metric = metrics,
      csv_val = csv_vals[metrics],
      ref_val = rds_vals[metrics],
      abs_diff = abs(csv_vals - rds_vals),
      pass = abs(csv_vals - rds_vals) <= tol,
      check.names = FALSE
    )
  }))
}

# ----------------------
# File Paths (relative to testthat.R)
# ----------------------
script_dir <- dirname(normalizePath(sys.frame(1)$ofile %||% "."))
project_root <- normalizePath(file.path(script_dir, ".."))

csv_dir <- file.path(project_root, "output")
rds_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "output/validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

cat("Project root: ", project_root, "\n")

if (length(csv_files) != 1 || length(rds_files) != 1) {
  stop("Expect exactly one CSV and one RDS file")
}

csv_file <- csv_files[1]
rds_file <- rds_files[1]
out_file <- file.path(output_dir, "validation_results.csv")

# ----------------------
# Run validation
# ----------------------
cmp <- validate_metrics(csv_file, rds_file, metrics, tol = 1e-3)
write.csv(cmp, out_file, row.names = FALSE, quote = FALSE)

# ----------------------
# Tests: one test per sample per metric
# ----------------------
for (s in unique(cmp$sample)) {
  for (m in metrics) {
    test_that(paste("sample:", s, "metric:", m), {
      expect_true(cmp$pass[cmp$sample == s & cmp$metric == m])
    })
  }
}

cat("Validation results written to:", out_file, "\n")


