# Code obtained from https://github.com/WEHIGenomicsRnD/G000204_duplex/blob/main/code/efficiency_nanoseq_functions.R

# Efficiency metric calculations
calculate_singletons <- function(rbs) {
    total_reads <- sum(rbs$x, rbs$y)
    singletons <- sum(rbs$x == 1 & rbs$y == 0 | rbs$x == 0 & rbs$y == 1)
    frac_singletons <- singletons / total_reads
    return(frac_singletons)
}

calculate_family_stats <- function(rbs) {
    rbs$size <- rbs$x + rbs$y
    return(c(total_families = nrow(rbs),
             family_mean = mean(rbs$size),
             family_median = median(rbs$size),
             family_max = max(rbs$size),
             families_gt1 = sum(rbs$x > 1 | rbs$y > 1),
             single_families = sum(rbs$x == 1 & rbs$y == 0 | rbs$x == 0 & rbs$y == 1),
             paired_families = sum(rbs$x > 0 & rbs$y > 0),
             paired_and_gt1 = sum(rbs$x > 1 & rbs$y > 1)))
}

calculate_metrics <- function(rbs) {
  metrics <- data.frame(sample = names(rbs))
  metrics$frac_singletons <- lapply(rbs, calculate_singletons) %>% unlist()
  metrics$efficiency <- lapply(rbs, calculate_efficiency) %>% unlist()
  metrics$drop_out_rate <- lapply(rbs, calculate_missed_fraction) %>% unlist()
  
  metrics <- lapply(rbs, function(x) calculate_gc(x, genomeFile = genomeFile, genome_max = genome_max)) %>%
    data.frame %>% t() %>%
    cbind(metrics, .)
  
  metrics <- lapply(rbs, calculate_family_stats) %>%
    data.frame %>% t() %>%
    cbind(metrics, .)
  return(metrics)
}



# --- Metric grouping / selection ---------

.individual_metrics <- c("frac_singletons", "efficiency", "drop_out_rate")

.metric_groups <- list(
  gc = c("gc_single", "gc_both", "gc_deviation"),
  family = c(
    "total_families", "family_mean", "family_median", "family_max",
    "families_gt1", "single_families", "paired_families", "paired_and_gt1"
  )
)

# Resolve --metrics into:
# - groups: grouped metrics to compute (gc/family)
# - individual: other metrics to compute individually (efficiency, drop_out_rate, frac_singletons)
#
# Rules:
# - empty / NULL -> compute all avaliable metrics
# - token "gc" or "family" -> compute that whole group
# - token is a individual metric name -> compute only that metric
# - token is a metric inside gc/family -> compute the whole group 
resolve_metric_selection <- function(metrics_arg = NULL) {
  if (is.null(metrics_arg) || !nzchar(metrics_arg)) {
    return(list(groups = names(.metric_groups), individual = .individual_metrics))
  }
  
  tokens <- trimws(unlist(strsplit(metrics_arg, ",")))
  tokens <- tokens[nzchar(tokens)]
  
  groups <- character(0)
  individual  <- character(0)
  
  for (tok in tokens) {
    if (tok %in% names(.metric_groups)) {
      groups <- union(groups, tok)
      next
    }
    if (tok %in% .individual_metrics) {
      individual <- union(individual, tok)
      next
    }
    
    # if metric name inside a group, map to its group
    hit <- names(Filter(function(v) tok %in% v, .metric_groups))
    if (length(hit) > 0) {
      groups <- union(groups, hit)
      next
    }
    
    stop("Unknown metric/group in --metrics: ", tok,
         "\nValid groups: ", paste(names(.metric_groups), collapse = ", "),
         "\nInidividual metrics: ", paste(.individual_metrics, collapse = ", "),
         "\nGrouped metrics: ", paste(unique(unlist(.metric_groups)), collapse = ", "))
  }
  
  list(groups = groups, individual = individual)
}

# Compute selected metrics (returns 1-row data.frame)
calculate_metrics_selected <- function(rbs,
                                       groups = c("gc", "family"),
                                       individual = character(0)) {
  metrics <- list()
  
  if ("frac_singletons" %in% individual) metrics$frac_singletons <- calculate_singletons(rbs)
  if ("efficiency" %in% individual)      metrics$efficiency      <- calculate_efficiency(rbs)
  if ("drop_out_rate" %in% individual)   metrics$drop_out_rate   <- calculate_missed_fraction(rbs)
  
  if ("gc" %in% groups) {
    gc_stats <- calculate_gc(rbs, genomeFile = genomeFile, genome_max = genome_max)
    metrics <- c(metrics, as.list(gc_stats))
  }
  
  if ("family" %in% groups) {
    fam_stats <- calculate_family_stats(rbs)
    metrics <- c(metrics, as.list(fam_stats))
  }
  
  as.data.frame(metrics, check.names = FALSE)
}



calc_metrics_new_rbs <- function(rinfo_dir, pattern="\\.txt.gz", cores=8,
                                 metrics_arg = "", genomeFile = NULL, genome_max = NULL) {
  groups <- resolve_metric_groups(metrics_arg)
  
  metrics <-
    list.files(rinfo_dir, full.names = TRUE, recursive = TRUE, pattern = pattern) %>%
    mclapply(., fread, mc.cores = cores) %>%
    mclapply(., function(one_rbs) {
      if ("gc" %in% groups) {
        if (is.null(genomeFile) || is.null(genome_max)) {
          stop("GC requested but genomeFile/genome_max not provided to calc_metrics_new_rbs().")
        }
        local_genomeFile <- genomeFile
        local_genome_max <- genome_max
        genomeFile <- local_genomeFile
        genome_max <- local_genome_max
      }
      calculate_metrics_selected(one_rbs, groups = groups)
    }, mc.cores = cores)
  
  return(metrics)
}




# functions below are adapted from
# R/efficiency_nanoseq.R and perl/efficiency_nanoseq.pl
# from https://github.com/cancerit/NanoSeq

# from cancerit/NanoSeq documentation:
# "This is the number of duplex bases divided by the number of sequenced bases."
calculate_efficiency <- function(rbs) {
    bases_ok_rbs <- nrow(rbs[rbs$x > 1 & rbs$y > 1,]) * ((rlen - skips) * 2)
    total_reads <- sum(c(rbs$x, rbs$y))
    bases_sequenced <- total_reads * rlen * 2
    eff <- bases_ok_rbs / bases_sequenced
    return(eff)
}

# from cancerit/NanoSeq documentation:
# "This shows the fraction of read bundles missing one of the two
# original strands beyond what would be expected under random sampling
# (assuming a binomial process).
calculate_missed_fraction <- function(rbs) {
    rbs <- data.frame(rbs)
    rbs$size <- rbs$x + rbs$y
    rbs$size <- pmin(rbs$size, 10)
    total_missed <- 0
    for(size in c(4:10)) {
        exp_orphan <- (0.5 ** size) * 2
        total_this_size = nrow(rbs[which(rbs$size == size),])
        if(total_this_size > 0) {
            with_both_strands <- nrow(rbs[which(rbs$size == size & rbs$x > 0 & rbs$y > 0),])
            obs_orphan <- 1 - with_both_strands / total_this_size
            missed <- (obs_orphan - exp_orphan) * total_this_size
            total_missed <- total_missed + missed
        }
    }
    total_missed_fraction = total_missed / nrow(rbs[which(rbs$size >= 4),])
    return(total_missed_fraction)
}

# from cancerit/NanoSeq documentation:
# The GC content of RBs with both strands and with just one strand.
# I return the difference between the two values.

calculate_gc <- function(rbs, sample_n = 10000, max_gap = 100000, genomeFile, genome_max) {
  rbs <- data.frame(rbs)
  colnames(rbs)[5:6] <- c("plus", "minus")
  
  # remove any chroms not in the sizes vector
  rbs <- rbs[rbs$chrom %in% names(genome_max), ]

  # compute end and drop invalid ranges early
  rbs$end <- rbs$mpos + rlen - skips

  rbs <- rbs[!is.na(rbs$pos) & !is.na(rbs$end) & rbs$pos > 0 & rbs$end >= rbs$pos, ]

  # remove records with mate positions exceeding genome max
  if (length(genome_max) > 1) {
    for (i in seq_along(genome_max)) {
      chrom <- names(genome_max)[i]
      chrom_max <- genome_max[[i]]
      rbs <- rbs[!(rbs$chrom == chrom & rbs$end > chrom_max), ]
    }
  } else {
    rbs <- rbs[rbs$end <= genome_max, ]
  }

  # remove records with large distances between mates
  rbs <- rbs[(rbs$end - rbs$pos) < max_gap, ]

  # remove zero-records (redundant but safe)
  rbs <- rbs[rbs$pos != 0, ]

  # split RBs
  rbs_both <- rbs[which(rbs$minus + rbs$plus >= 4 &
                        rbs$minus >= 2 &
                        rbs$plus  >= 2), ]
  rbs_single <- rbs[which(rbs$minus + rbs$plus > 4 &
                          (rbs$minus == 0 | rbs$plus == 0)), ]


  # guard: not enough RBs
  if (nrow(rbs_both) == 0 || nrow(rbs_single) == 0) {
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }

  # sample
  rbs_both <- rbs_both[sample(seq_len(nrow(rbs_both)),
                               min(sample_n, nrow(rbs_both))), ]
  rbs_single <- rbs_single[sample(seq_len(nrow(rbs_single)),
                                   min(sample_n, nrow(rbs_single))), ]

  # extract sequences
  seqs_both <- GRanges(rbs_both$chrom,
                       IRanges(start = rbs_both$pos,
                               end   = rbs_both$end)) %>%
    scanFa(genomeFile, .) %>% as.vector()

  seqs_single <- GRanges(rbs_single$chrom,
                         IRanges(start = rbs_single$pos,
                                 end   = rbs_single$end)) %>%
    scanFa(genomeFile, .) %>% as.vector()

  # drop NA sequences
  seqs_both   <- seqs_both[!is.na(seqs_both)]
  seqs_single <- seqs_single[!is.na(seqs_single)]

  if (length(seqs_both) == 0 || length(seqs_single) == 0) {
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }

  seqs_both_collapsed   <- paste(seqs_both, collapse = "")
  seqs_single_collapsed <- paste(seqs_single, collapse = "")

  if (is.na(seqs_both_collapsed) || is.na(seqs_single_collapsed) ||
      nchar(seqs_both_collapsed) == 0 || nchar(seqs_single_collapsed) == 0) {
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }

  gc_both   <- s2c(seqs_both_collapsed)   %>% GC()
  gc_single <- s2c(seqs_single_collapsed) %>% GC()

  c(gc_single = gc_single,
    gc_both = gc_both,
    gc_deviation = abs(gc_single - gc_both))
}

