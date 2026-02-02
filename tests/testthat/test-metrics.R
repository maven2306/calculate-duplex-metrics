library(testthat)

source("../../R/calculate_nanoseq_functions.R")

# ------------------------------------------------------------------------------
# Small rbs table
# ------------------------------------------------------------------------------

rbs_basic <- data.frame(
  chrom = rep("2e914854fabb46b9_1", 9),
  pos   = c(0, 0, 0, 0, 0, 0, 0, 0, 0),
  mpos  = c(2, 2, 2, 2, 2, 3, 3, 3, 4),
  umi   = c("ATC-GCC", "CAT-GCC", "GAT-GAT",
    "GTT-CCC", "TTT-GGA", "AAT-CCG",
    "GTT-ACC", "NTT-NCC", "TGC-GTT"),
  x = c(1, 1, 7, 7, 4, 1, 3, 1, 1),
  y = c(3, 0, 2, 5, 9, 0, 0, 0, 0),
  stringsAsFactors = FALSE )

rbs_empty <- rbs_basic[0,]

rlen  <- 100
skips <- 0
genomeFile <- "../../ref/Escherichia_coli_ATCC_10798.fasta"
genome_max <- c(chr1 = 1000000)
# ------------------------------------------------------------------------------
# calculate_singletons
# ------------------------------------------------------------------------------

test_that("calculate_singletons returns fraction in [0,1]", {
  val <- calculate_singletons(rbs_basic)
  expect_type(val, "double")
  expect_true(val >= 0 && val <= 1)
  expect_equal(val,0.0888888888888889, tolerance = 1e-3)
})

test_that("calculate_singletons returns NA with empty input", {
  val <- calculate_singletons(rbs_empty)
  expect_true(is.na(val))
})

# ------------------------------------------------------------------------------
# calculate_family_stats
# ------------------------------------------------------------------------------

test_that("calculate_family_stats returns named vector with expected names", {
  stats <- calculate_family_stats(rbs_basic)
  expect_type(stats, "double")
  expect_named(stats, c("total_families","family_mean","family_median","family_max",
      "families_gt1","single_families","paired_families","paired_and_gt1"))
})

test_that("calculate_family_stats returns correct known values", {
  stats <- calculate_family_stats(rbs_basic)
  expect_equal(stats[["total_families"]], 9)
  expect_equal(stats[["family_mean"]], 5)
  expect_equal(stats[["family_max"]], 13)
  expect_equal(stats[["single_families"]], 4)
  expect_equal(stats[["paired_families"]], 4)
})

# ------------------------------------------------------------------------------
# calculate_efficiency
# ------------------------------------------------------------------------------

test_that("calculate_efficiency returns valid output", {
  eff <- calculate_efficiency(rbs_basic, rlen = rlen, skips = skips)
  expect_type(eff, "double")
  expect_true(is.finite(eff))
  expect_true(eff >= 0 && eff <= 1)
  expect_equal(eff, 0.0667, tolerance = 1e-3)
})

test_that("calculate_efficiency errors on invalid rlen/skips", {
  expect_error(calculate_efficiency(rbs_basic, rlen = -1, skips = 0))
  expect_error(calculate_efficiency(rbs_basic, rlen = 10, skips = 10))
})

test_that("calculate_efficiency returns NA for zero reads", {
  eff <- calculate_efficiency(rbs_empty, rlen = rlen, skips = skips)
  expect_true(is.na(eff))
})
# ------------------------------------------------------------------------------
# calculate_missed_fraction
# ------------------------------------------------------------------------------

test_that("calculate_missed_fraction returns numeric or NA", {
  val <- calculate_missed_fraction(rbs_basic)
  expect_type(val, "double")
  expect_equal(val,-0.033203125)
})

test_that("calculate_missed_fraction returns NA with empty input", {
  val <- calculate_missed_fraction(rbs_empty)
  expect_true(is.na(val))
})

# ------------------------------------------------------------------------------
# calculate_gc
# ------------------------------------------------------------------------------

test_that("calculate_gc returns NA metrics if reference genome missing", {
  gc <- calculate_gc(
    rbs_basic,
    rlen = rlen,
    skips = skips,
    genomeFile = NULL,
    genome_max = NULL
  )
  expect_named(gc, c("gc_single","gc_both","gc_deviation"))
  expect_true(all(is.na(gc)))
})


test_that("calculate_gc returns expected GC metrics when reference is provided", {
  gc <- calculate_gc(rbs_basic,rlen = rlen, skips = skips, genomeFile = genomeFile, genome_max = genome_max)
  expect_named(gc, c("gc_single", "gc_both", "gc_deviation"))
  expect_true(all(is.na(gc) | (gc >= 0 & gc <= 1)))
})

# ------------------------------------------------------------------------------
# resolve_metric_selection
# ------------------------------------------------------------------------------

test_that("resolve_metric_selection default selects all", {
  sel <- resolve_metric_selection()
  expect_true(all(c("gc","family") %in% sel$groups))
  expect_true(all(c("frac_singletons","efficiency","drop_out_rate") %in% sel$individual))
})

test_that("resolve_metric_selection errors on unknown metric", {
  expect_error(
    resolve_metric_selection("not_a_metric"),
    "Unknown metric/group"
  )
})

# ------------------------------------------------------------------------------
# calculate_metrics_selected
# ------------------------------------------------------------------------------

test_that("calculate_metrics_selected returns 1-row data.frame", {
  res <- calculate_metrics_selected(
    rbs_basic,
    groups = "family",
    individual = "frac_singletons",
    rlen = rlen,
    skips = skips )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true("frac_singletons" %in% colnames(res))
  expect_true("family_mean" %in% colnames(res))
})

test_that("calculate_metrics_selected errors if GC selected without refernce genome", {
  expect_error(
    calculate_metrics_selected(
      rbs_basic,
      groups = "gc",
      individual = character(0),
      rlen = rlen,
      skips = skips
    ),
    "GC metrics requested"
  )
})
