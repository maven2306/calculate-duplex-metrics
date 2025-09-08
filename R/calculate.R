process_data <- function(input, output, sample=NULL, rfunc_dir, rlen=151, skips=5) {
  if (is.null(sample) || !nzchar(sample)) {
    sample <- sub("\\.txt(\\.gz)?$", "", basename(input), ignore.case = TRUE)
  }
  script <- file.path("R","calc_duplex_metrics.R")
  if (!file.exists(script)) return(list(success = FALSE, error = paste0("Missing: ", script)))
  
  status <- system2("Rscript", c(script, input, output, sample, rfunc_dir, rlen, skips))
  if (status != 0) return(list(success = FALSE, error = paste0("calc_duplex_metrics.R exit status ", status)))
  
  list(success = TRUE)
}
