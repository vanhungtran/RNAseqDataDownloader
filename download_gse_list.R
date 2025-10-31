#!/usr/bin/env Rscript
# Script to download a list of GSE datasets
# Usage: Rscript download_gse_list.R

library(RNAseqDataDownloader)

# Function to get the default GSE list
get_default_gse_list <- function() {
  c("GSE102628", "GSE102641", "GSE102725", "GSE103489", "GSE104509", "GSE106087",
    "GSE106992", "GSE107361", "GSE107871", "GSE109182", "GSE109248", "GSE111053",
    "GSE111054", "GSE111055", "GSE11307", "GSE114286", "GSE114729", "GSE116486",
    "GSE117239", "GSE117405", "GSE11903", "GSE120721", "GSE120899", "GSE121212",
    "GSE123785", "GSE123786", "GSE123787", "GSE124700", "GSE124701", "GSE130588",
    "GSE133385", "GSE133477", "GSE13355", "GSE14905", "GSE16161", "GSE18686",
    "GSE18948", "GSE20264", "GSE24767", "GSE26866", "GSE26952", "GSE2737",
    "GSE27887", "GSE30355", "GSE30768", "GSE30999", "GSE31652", "GSE32407",
    "GSE32473", "GSE32620", "GSE32924", "GSE34248", "GSE36381", "GSE36387",
    "GSE36842", "GSE38039", "GSE40033", "GSE40263", "GSE41662", "GSE41663",
    "GSE41664", "GSE41745", "GSE41905", "GSE42305", "GSE42632", "GSE47598",
    "GSE47751", "GSE47944", "GSE47965", "GSE48586", "GSE50598", "GSE50614",
    "GSE50790", "GSE51440", "GSE52361", "GSE52471", "GSE53431", "GSE53552",
    "GSE54456", "GSE55201", "GSE5667", "GSE57225", "GSE57376", "GSE57383",
    "GSE57386", "GSE57405", "GSE58121", "GSE58558", "GSE58749", "GSE59294",
    "GSE60481", "GSE60709", "GSE60971", "GSE61281", "GSE62408", "GSE63079",
    "GSE63741", "GSE63979", "GSE63980", "GSE121212", "GSE141570", "GSE65832",
    "GSE6601", "GSE66511", "GSE6710", "GSE67785", "GSE67853", "GSE68923",
    "GSE68924", "GSE68939", "GSE69967", "GSE72246", "GSE74697", "GSE75343",
    "GSE75890", "GSE77719", "GSE78023", "GSE78097", "GSE79704", "GSE80047",
    "GSE80429", "GSE82140", "GSE83582", "GSE83645", "GSE85034", "GSE86451",
    "GSE89725", "GSE92472", "GSE93423", "GSE99802", "GSE224783", "GSE140684",
    "GSE15719", "GSE95759", "GSE67785", "GSE157194", "GSE182740", "GSE261704",
    "GSE283265")
}

# Main download function
download_all_gse <- function(gse_list = NULL, 
                             output_dir = "geo_batch_data",
                             data_type = "matrix",
                             continue_on_error = TRUE,
                             save_log = TRUE) {
  
  if (is.null(gse_list)) {
    gse_list <- get_default_gse_list()
  }
  
  # Remove duplicates
  gse_list <- unique(gse_list)
  
  cat("========================================\n")
  cat("Batch GSE Download Starting\n")
  cat("Total datasets:", length(gse_list), "\n")
  cat("Output directory:", output_dir, "\n")
  cat("========================================\n\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize tracking
  results <- list(
    success = character(),
    failed = character(),
    skipped = character(),
    errors = list()
  )
  
  # Log file
  log_file <- NULL
  if (save_log) {
    log_file <- file.path(output_dir, paste0("download_log_", 
                          format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    cat("Log file:", log_file, "\n\n")
    writeLines(paste("Batch download started at:", Sys.time()), log_file)
    write(paste("Total datasets:", length(gse_list)), log_file, append = TRUE)
    write("", log_file, append = TRUE)
  }
  
  # Download each GSE
  for (i in seq_along(gse_list)) {
    gse_id <- gse_list[i]
    
    cat(sprintf("\n[%d/%d] Processing %s...\n", i, length(gse_list), gse_id))
    
    if (!is.null(log_file)) {
      write(sprintf("[%d/%d] %s - Started at %s", 
                   i, length(gse_list), gse_id, Sys.time()), 
           log_file, append = TRUE)
    }
    
    tryCatch({
      result <- download_geo(
        accession = gse_id,
        output_dir = output_dir,
        data_type = data_type,
        verbose = FALSE
      )
      
      if (result$success) {
        results$success <- c(results$success, gse_id)
        cat(sprintf("✓ SUCCESS: %s downloaded to %s\n", gse_id, result$data_path))
        
        if (!is.null(log_file)) {
          write(sprintf("  SUCCESS - Files: %d, Path: %s", 
                       length(result$files), result$data_path), 
               log_file, append = TRUE)
        }
      } else {
        results$failed <- c(results$failed, gse_id)
        cat(sprintf("✗ FAILED: %s\n", gse_id))
        
        if (!is.null(log_file)) {
          write("  FAILED", log_file, append = TRUE)
        }
      }
      
    }, error = function(e) {
      results$failed <<- c(results$failed, gse_id)
      results$errors[[gse_id]] <<- e$message
      
      cat(sprintf("✗ ERROR: %s - %s\n", gse_id, e$message))
      
      if (!is.null(log_file)) {
        write(sprintf("  ERROR: %s", e$message), log_file, append = TRUE)
      }
      
      if (!continue_on_error) {
        stop(sprintf("Download stopped at %s due to error: %s", gse_id, e$message))
      }
    })
    
    # Small delay to avoid overwhelming the server
    Sys.sleep(2)
  }
  
  # Summary
  cat("\n========================================\n")
  cat("Batch Download Complete\n")
  cat("========================================\n")
  cat(sprintf("Total processed: %d\n", length(gse_list)))
  cat(sprintf("Successful: %d\n", length(results$success)))
  cat(sprintf("Failed: %d\n", length(results$failed)))
  
  if (length(results$failed) > 0) {
    cat("\nFailed datasets:\n")
    cat(paste("-", results$failed, collapse = "\n"), "\n")
  }
  
  if (!is.null(log_file)) {
    write("", log_file, append = TRUE)
    write("========================================", log_file, append = TRUE)
    write(sprintf("Batch download completed at: %s", Sys.time()), log_file, append = TRUE)
    write(sprintf("Total processed: %d", length(gse_list)), log_file, append = TRUE)
    write(sprintf("Successful: %d", length(results$success)), log_file, append = TRUE)
    write(sprintf("Failed: %d", length(results$failed)), log_file, append = TRUE)
    
    if (length(results$failed) > 0) {
      write("\nFailed datasets:", log_file, append = TRUE)
      write(paste("-", results$failed), log_file, append = TRUE)
    }
    
    cat("\nLog saved to:", log_file, "\n")
  }
  
  # Return results invisibly
  invisible(results)
}

# Run if script is executed directly
if (!interactive()) {
  # Parse command line arguments if needed
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0) {
    output_dir <- args[1]
  } else {
    output_dir <- "geo_batch_data"
  }
  
  # Run the batch download
  results <- download_all_gse(
    gse_list = get_default_gse_list(),
    output_dir = output_dir,
    data_type = "matrix",
    continue_on_error = TRUE,
    save_log = TRUE
  )
}
