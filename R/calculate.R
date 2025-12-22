#process_data <- function(input, output, sample=NULL, rfunc_dir, rlen=151, skips=5) {
#  if (is.null(sample) || !nzchar(sample)) {
#    sample <- sub("\\.txt(\\.gz)?$", "", basename(input), ignore.case = TRUE)
#  }
#  script <- file.path("R","calc_duplex_metrics.R")
#  if (!file.exists(script)) return(list(success = FALSE, error = paste0("Missing: ", script)))
  
#  status <- system2("Rscript", c(script, input, output, sample, rfunc_dir, rlen, skips))
#  if (status != 0) return(list(success = FALSE, error = paste0("calc_duplex_metrics.R exit status ", status)))
  
#  list(success = TRUE)
#}

source("R/calc_duplex_metrics.R")

process_data <- function(
    input, output,
    sample = NULL,
    rlen = 151,
    skips = 5,
    ref_fasta = "",
    skip_gc = TRUE,
    metrics = "",
    cores = 1
) {
  # input can be a single file or a vector of files
  if (is.null(input) || length(input) == 0) {
    return(list(success = FALSE, error = "No input files provided"))
  }
  
  # validate cores
  if (is.na(cores) || cores < 1) {
    return(list(success = FALSE, error = "--cores must be >= 1"))
  }
  
  # ensure output dir exists
  odir <- dirname(output)
  if (!dir.exists(odir)) dir.create(odir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- single input file path ----
  if (length(input) == 1) {
    # defaults
    if (is.null(sample) || !nzchar(sample)) {
      sample <- sub("\\.txt(\\.gz)?$", "", basename(input), ignore.case = TRUE)
    }
    
    # validate input
    if (!file.exists(input)) {
      return(list(success = FALSE, error = paste0("Input not found: ", input)))
    }
    
    ok <- tryCatch({
      out_df <- calc_duplex_metrics_one_file_df(
        input       = input,
        sample      = sample,
        rlen        = rlen,
        skips       = skips,
        ref_fasta   = ref_fasta,
        skip_gc     = skip_gc,
        metrics_arg = metrics
      )
      write.csv(out_df, output, row.names = FALSE, quote = FALSE)
      TRUE
    }, error = function(e) e)
    
    if (inherits(ok, "error")) {
      return(list(success = FALSE, error = ok$message))
    }
    
    return(list(success = TRUE))
  }
  
  # ---- multiple input files ----
  # sample must be null/empty for multi-file mode 
  if (!is.null(sample) && nzchar(sample)) {
    return(list(success = FALSE, error = "--sample can only be used with a single input file"))
  }
  
  # validate all input files exist 
  missing <- input[!file.exists(input)]
  if (length(missing) > 0) {
    return(list(success = FALSE,
                error = paste0("Input file(s) not found:\n", paste(missing, collapse = "\n"))))
  }
  
  # compute combined long-format df, then write output once
  out_df <- tryCatch({
    calc_duplex_metrics_many_files_df(
      inputs     = input,
      rlen       = rlen,
      skips      = skips,
      ref_fasta  = ref_fasta,
      skip_gc    = skip_gc,
      metrics_arg = metrics,
      cores      = cores
    )
  }, error = function(e) e)
  
  if (inherits(out_df, "error")) {
    return(list(success = FALSE, error = out_df$message))
  }
  
  ok <- tryCatch({
    write.csv(out_df, output, row.names = FALSE, quote = FALSE)
    TRUE
  }, error = function(e) e)
  
  if (inherits(ok, "error")) {
    return(list(success = FALSE, error = paste0("Write failed: ", ok$message)))
  }
  
  list(success = TRUE)
}



