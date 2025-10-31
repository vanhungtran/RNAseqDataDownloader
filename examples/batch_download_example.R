# Example: Batch Download Multiple GSE Datasets
# This script demonstrates how to download multiple GSE datasets at once

library(RNAseqDataDownloader)

# Load the batch download script
source("download_gse_list.R")

# Example 1: Download all datasets from the default list
# ------------------------------------------------------
results <- download_all_gse(
  gse_list = NULL,  # Use default list
  output_dir = "geo_batch_data",
  data_type = "matrix",
  continue_on_error = TRUE,
  save_log = TRUE
)

# Example 2: Download a custom subset of datasets
# ------------------------------------------------------
custom_list <- c("GSE102628", "GSE102641", "GSE102725", "GSE103489", "GSE104509")

results <- download_all_gse(
  gse_list = custom_list,
  output_dir = "geo_custom_data",
  data_type = "matrix",
  continue_on_error = TRUE,
  save_log = TRUE
)

# Example 3: Download with supplementary files
# ------------------------------------------------------
# For datasets where you need the raw data files
download_with_supp <- function(gse_list, output_dir = "geo_with_supp") {
  
  results <- list(success = character(), failed = character())
  
  for (gse_id in gse_list) {
    cat("\nDownloading", gse_id, "with supplementary files...\n")
    
    tryCatch({
      result <- download_geo(
        accession = gse_id,
        output_dir = output_dir,
        data_type = "both",  # Get both matrix and supplementary files
        get_supplementary = TRUE,
        verbose = TRUE
      )
      
      if (result$success) {
        results$success <- c(results$success, gse_id)
      } else {
        results$failed <- c(results$failed, gse_id)
      }
    }, error = function(e) {
      cat("Error:", e$message, "\n")
      results$failed <<- c(results$failed, gse_id)
    })
    
    Sys.sleep(3)  # Longer delay for supplementary files
  }
  
  return(results)
}

# Use it:
# small_list <- c("GSE102628", "GSE102641")
# results <- download_with_supp(small_list)

# Example 4: Check which datasets downloaded successfully
# ------------------------------------------------------
check_downloads <- function(gse_list, base_dir = "geo_batch_data") {
  
  status <- data.frame(
    GSE_ID = character(),
    Downloaded = logical(),
    Path = character(),
    Files = integer(),
    stringsAsFactors = FALSE
  )
  
  for (gse_id in gse_list) {
    gse_dir <- file.path(base_dir, gse_id)
    
    if (dir.exists(gse_dir)) {
      files <- list.files(gse_dir, recursive = TRUE)
      status <- rbind(status, data.frame(
        GSE_ID = gse_id,
        Downloaded = TRUE,
        Path = gse_dir,
        Files = length(files),
        stringsAsFactors = FALSE
      ))
    } else {
      status <- rbind(status, data.frame(
        GSE_ID = gse_id,
        Downloaded = FALSE,
        Path = NA,
        Files = 0,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(status)
}

# Use it:
# status <- check_downloads(get_default_gse_list())
# print(status)
# cat("\nDownloaded:", sum(status$Downloaded), "/", nrow(status), "\n")

# Example 5: Resume failed downloads
# ------------------------------------------------------
resume_failed <- function(results, output_dir = "geo_batch_data") {
  
  if (length(results$failed) == 0) {
    cat("No failed downloads to resume.\n")
    return(NULL)
  }
  
  cat("Resuming", length(results$failed), "failed downloads...\n")
  
  return(download_all_gse(
    gse_list = results$failed,
    output_dir = output_dir,
    data_type = "matrix",
    continue_on_error = TRUE,
    save_log = TRUE
  ))
}

# Use it:
# results2 <- resume_failed(results)

# Example 6: Parallel download (use with caution!)
# ------------------------------------------------------
# library(parallel)
# 
# parallel_download <- function(gse_list, output_dir = "geo_parallel", n_cores = 2) {
#   
#   cl <- makeCluster(n_cores)
#   clusterEvalQ(cl, library(RNAseqDataDownloader))
#   
#   results <- parLapply(cl, gse_list, function(gse_id) {
#     tryCatch({
#       result <- download_geo(
#         accession = gse_id,
#         output_dir = output_dir,
#         data_type = "matrix",
#         verbose = FALSE
#       )
#       list(gse_id = gse_id, success = result$success, error = NULL)
#     }, error = function(e) {
#       list(gse_id = gse_id, success = FALSE, error = e$message)
#     })
#   })
#   
#   stopCluster(cl)
#   return(results)
# }

cat("\n=== Batch Download Functions Loaded ===\n")
cat("Available functions:\n")
cat("  - download_all_gse(gse_list, output_dir, ...)\n")
cat("  - download_with_supp(gse_list, output_dir)\n")
cat("  - check_downloads(gse_list, base_dir)\n")
cat("  - resume_failed(results, output_dir)\n")
cat("  - get_default_gse_list()\n")
