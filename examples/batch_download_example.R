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

# Example 7: Load downloaded matrices into a list
# ------------------------------------------------------
load_geo_batch_matrices <- function(base_dir = "geo_batch_data",
                                    pattern = "_series_matrix.txt.gz",
                                    simplify = TRUE,
                                    verbose = TRUE) {
  if (!dir.exists(base_dir)) {
    stop("Directory does not exist: ", base_dir)
  }

  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Package 'GEOquery' is required. Install it with install.packages('GEOquery').")
  }

  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package 'Biobase' is required. Install it with BiocManager::install('Biobase').")
  }

  matrix_files <- list.files(base_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)

  if (length(matrix_files) == 0) {
    if (verbose) {
      cat("No series matrix files found in", base_dir, "matching pattern", pattern, "\n")
    }
    return(list())
  }

  matrices <- vector("list", length(matrix_files))
  file_names <- character(length(matrix_files))

  for (i in seq_along(matrix_files)) {
    file_path <- matrix_files[[i]]
    base_name <- basename(file_path)
    cleaned_name <- tools::file_path_sans_ext(tools::file_path_sans_ext(base_name))
    file_names[[i]] <- cleaned_name

    if (verbose) {
      cat("Loading", cleaned_name, "from", file_path, "\n")
    }

    matrices[[i]] <- tryCatch({
      geo_obj <- GEOquery::getGEO(filename = file_path, GSEMatrix = TRUE)
      eset <- if (inherits(geo_obj, "ExpressionSet")) geo_obj else geo_obj[[1]]
      expr <- Biobase::exprs(eset)
      if (simplify) {
        expr <- as.matrix(expr)
      }
      expr
    }, error = function(e) {
      warning("Failed to load ", file_path, ": ", e$message)
      NULL
    })
  }

  names(matrices) <- file_names
  valid_entries <- !vapply(matrices, is.null, logical(1))

  if (verbose) {
    cat("Loaded", sum(valid_entries), "of", length(matrices), "matrix file(s).\n")
  }

  matrices[valid_entries]
}

list_GEO_data <- load_geo_batch_matrices()

# Example 8: Combine matrices into long table with GEO code
# ------------------------------------------------------
combine_geo_matrices <- function(matrices,
                                 drop_na = TRUE,
                                 verbose = TRUE) {
  if (!length(matrices)) {
    if (verbose) {
      cat("Empty matrix list provided.\n")
    }
    return(data.frame())
  }

  combined_list <- vector("list", length(matrices))
  keep_idx <- 1L

  for (name in names(matrices)) {
    mat <- matrices[[name]]
    if (is.null(mat)) {
      next
    }

    geo_code <- sub("_.*$", "", name)

    if (verbose) {
      cat("Combining matrix", name, "with GEO code", geo_code, "\n")
    }

    if (is.null(dim(mat)) || any(dim(mat) == 0)) {
      if (verbose) {
        cat("  Skipping", name, "(empty matrix).\n")
      }
      next
    }

    mat <- as.matrix(mat)
    n_features <- nrow(mat)
    n_samples <- ncol(mat)

    if (n_features == 0 || n_samples == 0) {
      if (verbose) {
        cat("  Skipping", name, "(no rows or columns).\n")
      }
      next
    }

    feature_names <- rownames(mat)
    if (is.null(feature_names) || all(feature_names == "")) {
      feature_names <- paste0("feature_", seq_len(n_features))
    }

    sample_names <- colnames(mat)
    if (is.null(sample_names) || all(sample_names == "")) {
      sample_names <- paste0("sample_", seq_len(n_samples))
    }

    df <- data.frame(
      Feature = rep(feature_names, times = n_samples),
      Sample = rep(sample_names, each = n_features),
      Value = as.vector(mat),
      GEOcode = geo_code,
      stringsAsFactors = FALSE
    )

    if (drop_na) {
      df <- df[!is.na(df$Value), , drop = FALSE]
    }

    if (!nrow(df)) {
      if (verbose) {
        cat("  Skipping", name, "(all values NA).\n")
      }
      next
    }

    combined_list[[keep_idx]] <- df
    keep_idx <- keep_idx + 1L
  }

  combined_list <- combined_list[seq_len(keep_idx - 1L)]

  if (!length(combined_list)) {
    if (verbose) {
      cat("No matrices contained data after processing.\n")
    }
    return(data.frame())
  }

  combined_df <- do.call(rbind, combined_list)
  rownames(combined_df) <- NULL
  combined_df <- combined_df[, c("GEOcode", "Feature", "Sample", "Value")]

  if (verbose) {
    cat("Combined rows:", nrow(combined_df), "from", length(combined_list), "matrix/matrices.\n")
  }

  combined_df
}

combined_GEO_data <- combine_geo_matrices(list_GEO_data)

cat("\n=== Batch Download Functions Loaded ===\n")
cat("Available functions:\n")
cat("  - download_all_gse(gse_list, output_dir, ...)\n")
cat("  - download_with_supp(gse_list, output_dir)\n")
cat("  - check_downloads(gse_list, base_dir)\n")
cat("  - resume_failed(results, output_dir)\n")
cat("  - get_default_gse_list()\n")
cat("  - load_geo_batch_matrices(base_dir, pattern, ...)\n")
cat("  - combine_geo_matrices(matrices, ...)\n")
