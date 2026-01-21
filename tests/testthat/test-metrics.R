library(testthat)
library(data.table)
library(GenomicRanges)
library(Biostrings)

source("../../R/efficiency_nanoseq_functions.R")

# Load the test RBS file
rbs <- fread("../../data/test.rinfo")
rbs_list <- list(sample1 = rbs)

# Define Global Variables Required by functions  
assign("genome_max",c(chr1 = 1000000), envir = .GlobalEnv)
assign("rlen", 100, envir = .GlobalEnv)
assign("skips", 0, envir = .GlobalEnv)


# Reference genome for GC tests
genomeFile <- "../../ref/Escherichia_coli_ATCC_10798.fasta"

# -----------------------------
# Unit tests
# -----------------------------

test_that("calculate_singletons returns a valid fraction", {
  frac <- calculate_singletons(rbs)
  expect_type(frac, "double")
  expect_true(is.finite(frac))
  expect_true(frac >= 0 && frac <= 1)
})

test_that("calculate_family_stats returns expected named vector", {
  stats <- calculate_family_stats(rbs)
  expect_type(stats, "double")
  expect_named(stats)
  expect_equal(unname(stats["total_families"]), nrow(rbs))
})

test_that("calculate_efficiency returns numeric efficiency", {
  eff <- calculate_efficiency(rbs)
  expect_type(eff, "double")
  expect_true(is.finite(eff))
  expect_true(eff >= 0 && eff <= 1)
})

test_that("calculate_missed_fraction returns numeric fraction", {
  miss <- calculate_missed_fraction(rbs)
  expect_type(miss, "double")
  expect_true(is.finite(miss))
})

test_that("calculate_gc returns expected GC metrics", {
  gc <- calculate_gc(rbs, genomeFile = genomeFile, genome_max = genome_max)
  expect_named(gc, c("gc_single", "gc_both", "gc_deviation"))
  expect_true(all(is.na(gc) | (gc >= 0 & gc <= 1)))
})

test_that("resolve_metric_selection handles empty input", {
  sel <- resolve_metric_selection()
  expect_named(sel, c("groups", "individual"))
  expect_true(length(sel$groups) > 0)
})

test_that("resolve_metric_selection errors on unknown metric", {
  expect_error(
    resolve_metric_selection("not_a_metric"),
    "Unknown metric"
  )
})

test_that("calculate_metrics_selected returns expected columns", {
  out <- calculate_metrics_selected(rbs,
                                    groups     = "family",
                                    individual = "frac_singletons")
  expect_s3_class(out, "data.frame")
  expect_true("frac_singletons" %in% colnames(out))
  expect_true("total_families" %in% colnames(out))
})

test_that("calculate_metrics returns a data.frame with expected metrics", {
  metrics <- calculate_metrics(rbs_list)
  expect_s3_class(metrics, "data.frame")

  expected_cols <- c(
    "frac_singletons",
    "efficiency",
    "drop_out_rate",
    "gc_single",
    "gc_both",
    "gc_deviation",
    "total_families",
    "family_mean",
    "family_median",
    "family_max"
  )

  expect_true(all(expected_cols %in% colnames(metrics)))
})
