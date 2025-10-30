#' Download data from GEO
#' 
#' Advanced function to download various types of GEO data including GSE (series),
#' GPL (platform), GSM (sample), and GDS (dataset) accessions. Supports multiple
#' data formats, RNA-seq quantifications, and includes comprehensive error handling.
#' Updated with latest GEOquery best practices from GitHub (2024-2025).
#' 
#' @param accession GEO accession number (e.g., "GSE12345", "GPL570", "GSM123456", "GDS1234")
#' @param output_dir Directory to save downloaded data (default: "geo_data")
#' @param data_type Type of data to download:
#'   \itemize{
#'     \item "matrix" - Expression matrix from series matrix files (GSE only, fastest)
#'     \item "raw" - Raw supplementary files (GSE/GSM)
#'     \item "full" - Complete GEO record with all metadata (SOFT format)
#'     \item "both" - Both matrix and raw data (GSE only)
#'     \item "rnaseq" - RNA-seq quantification data (if available from NCBI)
#'     \item "soft" - SOFT format files (legacy, slower)
#'   }
#' @param platform Character vector of GPL IDs to filter (for GSE with multiple platforms)
#' @param samples Character vector of GSM IDs to filter (subset of samples)
#' @param getGPL Logical, download platform annotation (default: TRUE)
#' @param AnnotGPL Logical, use annotation GPL with updated mappings (default: FALSE)
#' @param extract_archives Logical, automatically extract compressed archives (default: TRUE)
#' @param filter_files Regular expression to filter supplementary files
#' @param force_download Logical, re-download even if files exist (default: FALSE)
#' @param max_retries Integer, maximum download retry attempts (default: 3)
#' @param timeout Integer, download timeout in seconds (default: 300)
#' @param parseCharacteristics Logical, parse characteristics fields (default: TRUE)
#' @param verbose Logical, print detailed progress messages (default: TRUE)
#' @param ... Additional parameters passed to getGEO or getGEOSuppFiles
#' @return List containing:
#'   \itemize{
#'     \item data_path - Path to downloaded data
#'     \item geo_object - GEO object (ExpressionSet, SummarizedExperiment, or GEO class)
#'     \item metadata - Extracted metadata summary
#'     \item accession - Original accession number
#'     \item data_type - Type of data downloaded
#'     \item files - Vector of downloaded file paths
#'     \item platform_info - Platform information (if applicable)
#'     \item download_time - Timestamp of download
#'     \item success - Boolean indicating success
#'   }
#' @export
#' @examples
#' \dontrun{
#' # Download series matrix (fastest method)
#' result <- download_geo("GSE12345", data_type = "matrix")
#' 
#' # Download RNA-seq quantification data
#' result <- download_geo("GSE164073", data_type = "rnaseq")
#' 
#' # Download raw supplementary files with filtering
#' result <- download_geo("GSE12345", data_type = "raw", filter_files = "\\.txt\\.gz$")
#' 
#' # Download both matrix and supplementary files
#' result <- download_geo("GSE12345", data_type = "both")
#' 
#' # Download with annotation GPL (updated gene mappings)
#' result <- download_geo("GSE12345", AnnotGPL = TRUE)
#' 
#' # Download platform annotation
#' result <- download_geo("GPL570", data_type = "full")
#' 
#' # Download specific samples only
#' result <- download_geo("GSM123456", data_type = "full")
#' }
download_geo <- function(accession, 
                        output_dir = "geo_data",
                        data_type = c("matrix", "raw", "full", "both", "rnaseq", "soft"),
                        platform = NULL,
                        samples = NULL,
                        getGPL = TRUE,
                        AnnotGPL = FALSE,
                        extract_archives = TRUE,
                        filter_files = NULL,
                        force_download = FALSE,
                        max_retries = 3,
                        timeout = 300,
                        parseCharacteristics = TRUE,
                        verbose = TRUE,
                        ...) {
  
  # Validate inputs
  data_type <- match.arg(data_type)
  accession <- toupper(trimws(accession))
  
  if (!grepl("^(GSE|GPL|GSM|GDS)\\d+$", accession)) {
    stop("Invalid GEO accession format. Must be GSE, GPL, GSM, or GDS followed by numbers")
  }
  
  # Detect accession type
  acc_type <- substr(accession, 1, 3)
  
  # Set timeout option for downloads
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout))
  options(timeout = max(timeout, old_timeout))
  
  # Create output directory
  check_and_create_dir(output_dir)
  acc_dir <- file.path(output_dir, accession)
  check_and_create_dir(acc_dir)
  
  if (verbose) {
    message(paste0(strrep("=", 70)))
    message("Downloading GEO ", acc_type, " data: ", accession)
    message("Data type: ", data_type)
    message("Output directory: ", acc_dir)
    if (AnnotGPL) message("Using annotation GPL with updated gene mappings")
    message(paste0(strrep("=", 70)))
  }
  
  # Initialize result structure
  result <- list(
    accession = accession,
    data_type = data_type,
    data_path = acc_dir,
    geo_object = NULL,
    metadata = NULL,
    files = character(0),
    platform_info = NULL,
    download_time = Sys.time(),
    success = FALSE
  )
  
  # Download based on accession type and data type
  tryCatch({
    
    if (acc_type == "GSE") {
      result <- .download_gse(accession, acc_dir, data_type, platform, samples, 
                             getGPL, AnnotGPL, extract_archives, filter_files, 
                             force_download, max_retries, parseCharacteristics, 
                             verbose, ...)
    } else if (acc_type == "GPL") {
      result <- .download_gpl(accession, acc_dir, data_type, AnnotGPL, 
                             force_download, max_retries, verbose, ...)
    } else if (acc_type == "GSM") {
      result <- .download_gsm(accession, acc_dir, data_type, extract_archives,
                             force_download, max_retries, verbose, ...)
    } else if (acc_type == "GDS") {
      result <- .download_gds(accession, acc_dir, getGPL, force_download, 
                             max_retries, verbose, ...)
    }
    
    result$success <- TRUE
    
    if (verbose) {
      message("\n", strrep("=", 70))
      message("SUCCESS: Download completed for ", accession)
      message("Files saved to: ", result$data_path)
      message("Number of files: ", length(result$files))
      message(strrep("=", 70))
    }
    
    return(result)
    
  }, error = function(e) {
    warning("Failed to download GEO data for ", accession, ": ", e$message)
    result$error <- e$message
    return(result)
  })
}

#' Download GSE (Series) data
#' @keywords internal
.download_gse <- function(accession, acc_dir, data_type, platform, samples,
                          getGPL, AnnotGPL, extract_archives, filter_files,
                          force_download, max_retries, parseCharacteristics, 
                          verbose, ...) {
  
  result <- list(
    accession = accession,
    data_path = acc_dir,
    files = character(0),
    geo_object = NULL,
    metadata = NULL,
    platform_info = NULL
  )
  
  # Check for RNA-seq quantification data first
  if (data_type == "rnaseq") {
    if (verbose) message("\nChecking for RNA-seq quantification data...")
    
    attempt <- 0
    success <- FALSE
    
    while (attempt < max_retries && !success) {
      attempt <- attempt + 1
      
      tryCatch({
        # Try to get RNA-seq data using GEOquery's getRNASeqData
        if (requireNamespace("GEOquery", quietly = TRUE)) {
          se <- GEOquery::getRNASeqData(accession)
          result$geo_object <- se
          result$metadata <- list(
            type = "SummarizedExperiment",
            n_samples = ncol(se),
            n_features = nrow(se),
            assays = names(SummarizedExperiment::assays(se)),
            genome_info = S4Vectors::metadata(se)
          )
          
          # Save the SummarizedExperiment object
          se_file <- file.path(acc_dir, paste0(accession, "_rnaseq_se.rds"))
          saveRDS(se, se_file)
          result$files <- c(result$files, se_file)
          
          if (verbose) {
            message("RNA-seq quantification data downloaded successfully")
            message("Samples: ", ncol(se))
            message("Features: ", nrow(se))
          }
          
          success <- TRUE
        } else {
          stop("GEOquery package required for RNA-seq data download")
        }
        
      }, error = function(e) {
        if (attempt >= max_retries) {
          warning("No RNA-seq quantification data available for ", accession)
          warning("Error: ", e$message)
          if (verbose) {
            message("Falling back to standard matrix download...")
          }
          # Fall back to matrix download
          data_type <<- "matrix"
        } else {
          Sys.sleep(2 ^ attempt)
        }
      })
    }
    
    if (success) {
      return(result)
    }
  }
  
  # Download series matrix (GSEMatrix files - fastest method)
  if (data_type %in% c("matrix", "both", "full", "soft")) {
    if (verbose) message("\n[1/", ifelse(data_type == "both", "2", "1"), "] Downloading series matrix files...")
    
    attempt <- 0
    success <- FALSE
    
    while (attempt < max_retries && !success) {
      attempt <- attempt + 1
      if (verbose && attempt > 1) message("Retry attempt ", attempt, "/", max_retries)
      
      tryCatch({
        # Use GSEMatrix=TRUE for fast parsing (default in GEOquery 2.99.0+)
        use_gsematrix <- data_type != "soft"
        
        gse <- GEOquery::getGEO(
          accession, 
          destdir = acc_dir, 
          GSEMatrix = use_gsematrix,
          getGPL = getGPL,
          AnnotGPL = AnnotGPL,
          parseCharacteristics = parseCharacteristics,
          ...
        )
        
        # Handle list of ExpressionSets (multiple platforms)
        if (is.list(gse) && length(gse) > 0) {
          if (length(gse) > 1 && verbose) {
            message("Found ", length(gse), " platform(s) in this series:")
            for (i in seq_along(gse)) {
              plat <- Biobase::annotation(gse[[i]])
              n_samples <- ncol(gse[[i]])
              message("  [", i, "] ", plat, " (", n_samples, " samples)")
            }
          }
          
          # Filter by platform if specified
          if (!is.null(platform)) {
            platform_names <- sapply(gse, function(x) Biobase::annotation(x))
            matched_idx <- which(platform_names %in% platform)
            if (length(matched_idx) > 0) {
              gse <- gse[matched_idx]
              if (verbose) {
                message("Filtered to platform(s): ", paste(platform, collapse = ", "))
              }
            } else {
              warning("Specified platform(s) not found. Available: ", 
                     paste(platform_names, collapse = ", "))
            }
          }
        }
        
        result$geo_object <- gse
        result$metadata <- .extract_gse_metadata(gse, samples, verbose)
        
        # Get platform information
        if (is.list(gse) && length(gse) > 0) {
          result$platform_info <- lapply(gse, function(x) {
            list(
              platform_id = Biobase::annotation(x),
              n_features = nrow(x),
              feature_names = head(Biobase::featureNames(x), 10)
            )
          })
        }
        
        # List matrix files
        matrix_files <- list.files(
          acc_dir, 
          pattern = paste0(accession, ".*\\.txt(\\.gz)?$"), 
          full.names = TRUE
        )
        result$files <- c(result$files, matrix_files)
        
        success <- TRUE
        if (verbose) message("Series matrix downloaded successfully")
        
      }, error = function(e) {
        if (attempt >= max_retries) {
          stop("Failed to download series matrix after ", max_retries, 
               " attempts: ", e$message)
        }
        Sys.sleep(2 ^ attempt)  # Exponential backoff
      })
    }
  }
  
  # Download supplementary files
  if (data_type %in% c("raw", "both", "full")) {
    if (verbose) {
      step <- ifelse(data_type == "both", "[2/2]", "[1/1]")
      message("\n", step, " Downloading supplementary files...")
    }
    
    attempt <- 0
    success <- FALSE
    
    while (attempt < max_retries && !success) {
      attempt <- attempt + 1
      if (verbose && attempt > 1) message("Retry attempt ", attempt, "/", max_retries)
      
      tryCatch({
        # Use filter_regex parameter for file filtering
        supp_files <- GEOquery::getGEOSuppFiles(
          accession, 
          baseDir = dirname(acc_dir),
          makeDirectory = FALSE,
          fetch_files = TRUE,
          filter_regex = filter_files,
          ...
        )
        
        if (!is.null(supp_files) && nrow(supp_files) > 0) {
          # Extract archives if requested
          if (extract_archives) {
            .extract_geo_archives(rownames(supp_files), verbose)
          }
          
          result$files <- c(result$files, rownames(supp_files))
          if (verbose) {
            message("Downloaded ", nrow(supp_files), " supplementary file(s)")
          }
        } else {
          if (verbose) message("No supplementary files available")
        }
        
        success <- TRUE
        
      }, error = function(e) {
        if (attempt >= max_retries) {
          warning("Failed to download supplementary files after ", max_retries, 
                 " attempts: ", e$message)
          break
        } else {
          Sys.sleep(2 ^ attempt)
        }
      })
      if (attempt >= max_retries) success <- TRUE  # Continue even if suppl files fail
    }
  }
  
  return(result)
}

#' Download GPL (Platform) data
#' @keywords internal
.download_gpl <- function(accession, acc_dir, data_type, AnnotGPL, force_download, 
                         max_retries, verbose, ...) {
  
  if (verbose) {
    message("Downloading platform annotation: ", accession)
    if (AnnotGPL) {
      message("Using annotation GPL with updated Entrez Gene mappings")
    }
  }
  
  result <- list(
    accession = accession,
    data_path = acc_dir,
    files = character(0),
    geo_object = NULL,
    metadata = NULL,
    platform_info = NULL
  )
  
  attempt <- 0
  success <- FALSE
  
  while (attempt < max_retries && !success) {
    attempt <- attempt + 1
    
    tryCatch({
      gpl <- GEOquery::getGEO(accession, destdir = acc_dir, AnnotGPL = AnnotGPL, ...)
      result$geo_object <- gpl
      
      # Extract platform metadata
      meta <- GEOquery::Meta(gpl)
      table_data <- GEOquery::Table(gpl)
      
      result$metadata <- list(
        platform_id = accession,
        title = meta$title,
        organism = meta$organism,
        technology = meta$technology,
        manufacturer = meta$manufacturer,
        distribution = meta$distribution,
        n_probes = nrow(table_data),
        columns = colnames(table_data),
        submission_date = meta$submission_date,
        last_update = meta$last_update_date
      )
      
      result$platform_info <- result$metadata
      
      # Save annotation table
      annotation_file <- file.path(acc_dir, paste0(accession, "_annotation.txt"))
      write.table(table_data, annotation_file, sep = "\t", 
                 row.names = FALSE, quote = FALSE)
      result$files <- c(result$files, annotation_file)
      
      # Also save as RDS for easy loading
      rds_file <- file.path(acc_dir, paste0(accession, "_platform.rds"))
      saveRDS(gpl, rds_file)
      result$files <- c(result$files, rds_file)
      
      if (verbose) {
        message("Platform: ", result$metadata$title)
        message("Organism: ", result$metadata$organism)
        message("Technology: ", result$metadata$technology)
        message("Probes/Features: ", result$metadata$n_probes)
        if (AnnotGPL) {
          message("Annotation includes updated gene mappings")
        }
      }
      
      success <- TRUE
      
    }, error = function(e) {
      if (attempt >= max_retries) {
        stop("Failed to download platform data: ", e$message)
      }
      Sys.sleep(2 ^ attempt)
    })
  }
  
  return(result)
}

#' Download GSM (Sample) data
#' @keywords internal
.download_gsm <- function(accession, acc_dir, data_type, extract_archives,
                         force_download, max_retries, verbose, ...) {
  
  if (verbose) message("Downloading sample data: ", accession)
  
  result <- list(
    accession = accession,
    data_path = acc_dir,
    files = character(0),
    geo_object = NULL,
    metadata = NULL
  )
  
  attempt <- 0
  success <- FALSE
  
  while (attempt < max_retries && !success) {
    attempt <- attempt + 1
    
    tryCatch({
      # Download GSM metadata
      gsm <- GEOquery::getGEO(accession, destdir = acc_dir, ...)
      result$geo_object <- gsm
      
      # Extract sample metadata
      result$metadata <- list(
        sample_id = accession,
        title = GEOquery::Meta(gsm)$title,
        source = GEOquery::Meta(gsm)$source_name_ch1,
        organism = GEOquery::Meta(gsm)$organism_ch1,
        platform = GEOquery::Meta(gsm)$platform_id,
        series = GEOquery::Meta(gsm)$series_id,
        supplementary_files = GEOquery::Meta(gsm)$supplementary_file
      )
      
      if (verbose) {
        message("Sample: ", result$metadata$title)
        message("Series: ", result$metadata$series)
        message("Platform: ", result$metadata$platform)
      }
      
      # Download supplementary files if requested
      if (data_type %in% c("raw", "full")) {
        supp_files <- GEOquery::getGEOSuppFiles(accession, baseDir = dirname(acc_dir), ...)
        if (nrow(supp_files) > 0) {
          result$files <- c(result$files, rownames(supp_files))
          
          if (extract_archives) {
            .extract_geo_archives(rownames(supp_files), verbose)
          }
        }
      }
      
      success <- TRUE
      
    }, error = function(e) {
      if (attempt >= max_retries) {
        stop("Failed to download sample data: ", e$message)
      }
      Sys.sleep(2 ^ attempt)
    })
  }
  
  return(result)
}

#' Download GDS (Dataset) data
#' @keywords internal
.download_gds <- function(accession, acc_dir, getGPL, force_download, max_retries, verbose, ...) {
  
  if (verbose) message("Downloading GEO dataset: ", accession)
  
  result <- list(
    accession = accession,
    data_path = acc_dir,
    files = character(0),
    geo_object = NULL,
    metadata = NULL
  )
  
  attempt <- 0
  success <- FALSE
  
  while (attempt < max_retries && !success) {
    attempt <- attempt + 1
    
    tryCatch({
      gds <- GEOquery::getGEO(accession, destdir = acc_dir, ...)
      result$geo_object <- gds
      
      # Extract GDS metadata
      meta <- GEOquery::Meta(gds)
      result$metadata <- list(
        dataset_id = accession,
        title = meta$title,
        description = meta$description,
        platform = meta$platform,
        type = meta$type,
        n_samples = ncol(GEOquery::Table(gds)) - 2,
        update_date = meta$update_date
      )
      
      # Convert to ExpressionSet
      eset <- GEOquery::GDS2eSet(gds, do.log2 = FALSE, getGPL = getGPL)
      
      # Save expression data
      expr_file <- file.path(acc_dir, paste0(accession, "_expression.txt"))
      write.table(Biobase::exprs(eset), expr_file, sep = "\t", quote = FALSE)
      result$files <- c(result$files, expr_file)
      
      # Save phenotype data
      pheno_file <- file.path(acc_dir, paste0(accession, "_phenotype.txt"))
      write.table(Biobase::pData(eset), pheno_file, sep = "\t", quote = FALSE)
      result$files <- c(result$files, pheno_file)
      
      # Save as RDS
      rds_file <- file.path(acc_dir, paste0(accession, "_eset.rds"))
      saveRDS(eset, rds_file)
      result$files <- c(result$files, rds_file)
      
      if (verbose) {
        message("Dataset: ", result$metadata$title)
        message("Platform: ", result$metadata$platform)
        message("Samples: ", result$metadata$n_samples)
      }
      
      success <- TRUE
      
    }, error = function(e) {
      if (attempt >= max_retries) {
        stop("Failed to download GDS data: ", e$message)
      }
      Sys.sleep(2 ^ attempt)
    })
  }
  
  return(result)
}

#' Extract archives (tar.gz, gz, zip)
#' @keywords internal
.extract_geo_archives <- function(file_paths, verbose = TRUE) {
  for (file_path in file_paths) {
    if (!file.exists(file_path)) next
    
    if (grepl("\\.tar\\.gz$|\\.tgz$", file_path)) {
      if (verbose) message("Extracting: ", basename(file_path))
      utils::untar(file_path, exdir = dirname(file_path))
    } else if (grepl("\\.gz$", file_path) && !grepl("\\.tar\\.gz$", file_path)) {
      if (verbose) message("Extracting: ", basename(file_path))
      R.utils::gunzip(file_path, remove = FALSE, overwrite = TRUE)
    } else if (grepl("\\.zip$", file_path)) {
      if (verbose) message("Extracting: ", basename(file_path))
      utils::unzip(file_path, exdir = dirname(file_path))
    } else if (grepl("\\.tar$", file_path)) {
      if (verbose) message("Extracting: ", basename(file_path))
      utils::untar(file_path, exdir = dirname(file_path))
    }
  }
}

#' Extract comprehensive metadata from GSE
#' @keywords internal
.extract_gse_metadata <- function(gse_object, samples = NULL, verbose = TRUE) {
  
  metadata <- list()
  
  # Handle list of ExpressionSets (multiple platforms)
  if (is.list(gse_object)) {
    metadata$n_platforms <- length(gse_object)
    metadata$platforms <- list()
    
    for (i in seq_along(gse_object)) {
      eset <- gse_object[[i]]
      platform_id <- Biobase::annotation(eset)
      
      pheno <- Biobase::pData(eset)
      
      # Filter samples if specified
      if (!is.null(samples)) {
        sample_ids <- rownames(pheno)
        matched_samples <- intersect(samples, sample_ids)
        if (length(matched_samples) > 0) {
          pheno <- pheno[matched_samples, , drop = FALSE]
          if (verbose) {
            message("Platform ", platform_id, ": Filtered to ", 
                   length(matched_samples), " sample(s)")
          }
        }
      }
      
      metadata$platforms[[platform_id]] <- list(
        platform_id = platform_id,
        n_samples = nrow(pheno),
        n_features = nrow(Biobase::exprs(eset)),
        sample_ids = rownames(pheno),
        phenotype_variables = colnames(pheno)
      )
    }
  } else {
    # Single ExpressionSet
    pheno <- Biobase::pData(gse_object)
    
    if (!is.null(samples)) {
      sample_ids <- rownames(pheno)
      matched_samples <- intersect(samples, sample_ids)
      if (length(matched_samples) > 0) {
        pheno <- pheno[matched_samples, , drop = FALSE]
      }
    }
    
    metadata$platform_id <- Biobase::annotation(gse_object)
    metadata$n_samples <- nrow(pheno)
    metadata$n_features <- nrow(Biobase::exprs(gse_object))
    metadata$sample_ids <- rownames(pheno)
    metadata$phenotype_variables <- colnames(pheno)
  }
  
  return(metadata)
}



#' Extract comprehensive metadata from GEO object
#' 
#' Advanced metadata extraction supporting GSE, GPL, GSM, and GDS objects.
#' Extracts phenotype data, platform information, experimental design,
#' and creates summary statistics.
#' 
#' @param geo_object GEO object from getGEO (ExpressionSet, list of ExpressionSets, 
#'   GPL, GSM, or GDS object)
#' @param extract_pheno Logical, extract full phenotype data (default: TRUE)
#' @param extract_features Logical, extract feature/probe information (default: FALSE)
#' @param simplify Logical, simplify column names and remove redundant info (default: TRUE)
#' @param save_to File path to save metadata as TSV (optional)
#' @return Data frame or list with comprehensive metadata including:
#'   \itemize{
#'     \item phenotype_data - Sample phenotype/clinical data
#'     \item platform_info - Platform annotation details
#'     \item series_info - Series-level metadata
#'     \item summary - Summary statistics
#'     \item expression_data - Expression matrix (if requested)
#'   }
#' @export
#' @examples
#' \dontrun{
#' # Download and extract metadata
#' gse <- GEOquery::getGEO("GSE12345")
#' metadata <- extract_geo_metadata(gse)
#' 
#' # Access phenotype data
#' pheno <- metadata$phenotype_data
#' 
#' # Save metadata to file
#' metadata <- extract_geo_metadata(gse, save_to = "metadata.tsv")
#' }
extract_geo_metadata <- function(geo_object, 
                                extract_pheno = TRUE,
                                extract_features = FALSE,
                                simplify = TRUE,
                                save_to = NULL) {
  
  if (is.null(geo_object)) {
    stop("geo_object cannot be NULL")
  }
  
  # Detect object type
  obj_class <- class(geo_object)
  
  result <- list(
    object_type = obj_class,
    extraction_time = Sys.time(),
    phenotype_data = NULL,
    platform_info = NULL,
    series_info = NULL,
    expression_data = NULL,
    summary = list()
  )
  
  # Handle different GEO object types
  if (inherits(geo_object, "list")) {
    # List of ExpressionSets (multiple platforms)
    result <- .extract_gse_list_metadata(geo_object, extract_pheno, 
                                        extract_features, simplify)
  } else if (inherits(geo_object, "ExpressionSet")) {
    # Single ExpressionSet
    result <- .extract_expressionset_metadata(geo_object, extract_pheno, 
                                             extract_features, simplify)
  } else if (inherits(geo_object, "GPL")) {
    # Platform object
    result <- .extract_gpl_metadata(geo_object, extract_features)
  } else if (inherits(geo_object, "GSM")) {
    # Sample object
    result <- .extract_gsm_metadata(geo_object)
  } else if (inherits(geo_object, "GDS")) {
    # Dataset object
    result <- .extract_gds_metadata(geo_object, extract_pheno)
  } else {
    stop("Unsupported GEO object type: ", paste(obj_class, collapse = ", "),
         "\nSupported types: ExpressionSet, list of ExpressionSets, GPL, GSM, GDS")
  }
  
  # Save to file if requested
  if (!is.null(save_to) && !is.null(result$phenotype_data)) {
    write.table(result$phenotype_data, save_to, sep = "\t", 
               row.names = TRUE, quote = FALSE)
    message("Metadata saved to: ", save_to)
  }
  
  return(result)
}

#' Extract metadata from list of ExpressionSets
#' @keywords internal
.extract_gse_list_metadata <- function(geo_list, extract_pheno, extract_features, simplify) {
  
  n_platforms <- length(geo_list)
  
  result <- list(
    object_type = "GSE (multiple platforms)",
    n_platforms = n_platforms,
    platforms = list(),
    phenotype_data = list(),
    series_info = list(),
    summary = list()
  )
  
  # Extract series-level information from first platform
  if (n_platforms > 0) {
    eset <- geo_list[[1]]
    exp_data <- Biobase::experimentData(eset)
    
    result$series_info <- list(
      title = exp_data@title,
      abstract = exp_data@abstract,
      url = exp_data@url,
      pubmed_id = exp_data@pubMedIds
    )
  }
  
  # Process each platform
  all_pheno <- list()
  
  for (i in seq_along(geo_list)) {
    eset <- geo_list[[i]]
    platform_id <- Biobase::annotation(eset)
    
    message("Processing platform ", i, "/", n_platforms, ": ", platform_id)
    
    # Extract phenotype data
    if (extract_pheno) {
      pheno <- Biobase::pData(eset)
      
      if (simplify) {
        pheno <- .simplify_phenotype_data(pheno)
      }
      
      # Add platform identifier
      pheno$platform <- platform_id
      all_pheno[[platform_id]] <- pheno
    }
    
    # Platform-specific info
    result$platforms[[platform_id]] <- list(
      platform_id = platform_id,
      n_samples = ncol(eset),
      n_features = nrow(eset),
      sample_ids = colnames(eset),
      feature_ids = if(nrow(eset) > 0) rownames(eset)[seq_len(min(10, nrow(eset)))] else character(0),
      phenotype_variables = if(extract_pheno) colnames(Biobase::pData(eset)) else NULL
    )
    
    # Extract feature data if requested
    if (extract_features) {
      fdata <- Biobase::fData(eset)
      if (nrow(fdata) > 0) {
        result$platforms[[platform_id]]$feature_data <- head(fdata, 100)
      }
    }
  }
  
  # Combine phenotype data from all platforms
  if (extract_pheno && length(all_pheno) > 0) {
    result$phenotype_data <- do.call(rbind, all_pheno)
  }
  
  # Summary statistics
  result$summary <- list(
    total_platforms = n_platforms,
    total_samples = sum(sapply(geo_list, ncol)),
    total_features = sapply(geo_list, nrow),
    platform_ids = sapply(geo_list, Biobase::annotation)
  )
  
  return(result)
}

#' Extract metadata from single ExpressionSet
#' @keywords internal
.extract_expressionset_metadata <- function(eset, extract_pheno, extract_features, simplify) {
  
  result <- list(
    object_type = "ExpressionSet",
    platform_info = list(),
    phenotype_data = NULL,
    series_info = list(),
    summary = list()
  )
  
  # Series information
  exp_data <- Biobase::experimentData(eset)
  result$series_info <- list(
    title = exp_data@title,
    abstract = exp_data@abstract,
    url = exp_data@url,
    pubmed_id = exp_data@pubMedIds,
    lab = exp_data@lab,
    contact = exp_data@contact
  )
  
  # Platform information
  platform_id <- Biobase::annotation(eset)
  result$platform_info <- list(
    platform_id = platform_id,
    n_samples = ncol(eset),
    n_features = nrow(eset),
    sample_ids = colnames(eset),
    feature_ids = head(rownames(eset), 100)
  )
  
  # Phenotype data
  if (extract_pheno) {
    pheno <- Biobase::pData(eset)
    
    if (simplify) {
      pheno <- .simplify_phenotype_data(pheno)
    }
    
    result$phenotype_data <- pheno
    
    # Identify potential grouping variables
    result$summary$potential_groups <- .identify_grouping_variables(pheno)
  }
  
  # Feature data
  if (extract_features) {
    fdata <- Biobase::fData(eset)
    if (nrow(fdata) > 0) {
      result$platform_info$feature_annotation <- head(fdata, 100)
      result$summary$annotation_columns <- colnames(fdata)
    }
  }
  
  # Expression data summary
  expr <- Biobase::exprs(eset)
  result$summary$expression_stats <- list(
    min = min(expr, na.rm = TRUE),
    max = max(expr, na.rm = TRUE),
    mean = mean(expr, na.rm = TRUE),
    median = median(expr, na.rm = TRUE),
    na_count = sum(is.na(expr)),
    na_percentage = sum(is.na(expr)) / length(expr) * 100
  )
  
  return(result)
}

#' Extract metadata from GPL object
#' @keywords internal
.extract_gpl_metadata <- function(gpl_obj, extract_features) {
  
  meta <- GEOquery::Meta(gpl_obj)
  table_data <- GEOquery::Table(gpl_obj)
  
  result <- list(
    object_type = "GPL",
    platform_info = list(
      platform_id = meta$geo_accession,
      title = meta$title,
      organism = meta$organism,
      technology = meta$technology,
      distribution = meta$distribution,
      manufacturer = meta$manufacturer,
      manufacture_protocol = meta$manufacture_protocol,
      description = meta$description,
      n_probes = nrow(table_data),
      annotation_columns = colnames(table_data)
    ),
    summary = list(
      submission_date = meta$submission_date,
      last_update_date = meta$last_update_date,
      data_row_count = meta$data_row_count
    )
  )
  
  # Include feature annotation if requested
  if (extract_features && nrow(table_data) > 0) {
    result$feature_annotation <- table_data
  }
  
  return(result)
}

#' Extract metadata from GSM object
#' @keywords internal
.extract_gsm_metadata <- function(gsm_obj) {
  
  meta <- GEOquery::Meta(gsm_obj)
  
  result <- list(
    object_type = "GSM",
    sample_info = list(
      sample_id = meta$geo_accession,
      title = meta$title,
      source = meta$source_name_ch1,
      organism = meta$organism_ch1,
      molecule = meta$molecule_ch1,
      label = meta$label_ch1,
      platform = meta$platform_id,
      series = meta$series_id,
      description = meta$description,
      data_processing = meta$data_processing,
      supplementary_files = meta$supplementary_file
    ),
    characteristics = list(),
    summary = list(
      submission_date = meta$submission_date,
      last_update_date = meta$last_update_date
    )
  )
  
  # Extract characteristics
  char_cols <- grep("^characteristics_ch", names(meta), value = TRUE)
  if (length(char_cols) > 0) {
    result$characteristics <- meta[char_cols]
  }
  
  # Extract expression data if available
  if (!is.null(GEOquery::Table(gsm_obj))) {
    table_data <- GEOquery::Table(gsm_obj)
    result$expression_data <- table_data
    result$summary$n_features <- nrow(table_data)
  }
  
  return(result)
}

#' Extract metadata from GDS object
#' @keywords internal
.extract_gds_metadata <- function(gds_obj, extract_pheno) {
  
  meta <- GEOquery::Meta(gds_obj)
  table_data <- GEOquery::Table(gds_obj)
  
  result <- list(
    object_type = "GDS",
    dataset_info = list(
      dataset_id = meta$dataset_id,
      title = meta$title,
      description = meta$description,
      type = meta$type,
      platform = meta$platform,
      reference_series = meta$reference_series,
      n_samples = ncol(table_data) - 2,  # Exclude ID and IDENTIFIER columns
      n_features = nrow(table_data),
      sample_organism = meta$sample_organism
    ),
    summary = list(
      update_date = meta$update_date,
      pubmed_id = meta$pubmed_id,
      order = meta$order
    )
  )
  
  # Extract subset information
  if (!is.null(meta$subset_description)) {
    result$subsets <- list(
      descriptions = meta$subset_description,
      sample_ids = meta$subset_sample_id
    )
  }
  
  # Phenotype data from GDS
  if (extract_pheno) {
    # Convert GDS to ExpressionSet to get structured phenotype data
    tryCatch({
      eset <- GEOquery::GDS2eSet(gds_obj, do.log2 = FALSE)
      result$phenotype_data <- Biobase::pData(eset)
    }, error = function(e) {
      message("Could not extract phenotype data: ", e$message)
    })
  }
  
  return(result)
}

#' Simplify phenotype data by cleaning column names and values
#' @keywords internal
.simplify_phenotype_data <- function(pheno) {
  
  # Remove completely empty columns
  empty_cols <- sapply(pheno, function(x) all(is.na(x) | x == "" | x == "NA"))
  if (any(empty_cols)) {
    pheno <- pheno[, !empty_cols, drop = FALSE]
  }
  
  # Remove duplicate columns (same content, different names)
  if (ncol(pheno) > 1) {
    unique_cols <- !duplicated(t(pheno))
    pheno <- pheno[, unique_cols, drop = FALSE]
  }
  
  # Simplify column names
  colnames(pheno) <- gsub("^characteristics_ch1\\.", "", colnames(pheno))
  colnames(pheno) <- gsub("^characteristics_ch1$", "characteristics", colnames(pheno))
  colnames(pheno) <- gsub("[:_]+", "_", colnames(pheno))
  colnames(pheno) <- gsub("^_|_$", "", colnames(pheno))
  
  return(pheno)
}

#' Identify potential grouping variables in phenotype data
#' @keywords internal
.identify_grouping_variables <- function(pheno) {
  
  grouping_vars <- list()
  
  for (col in colnames(pheno)) {
    col_data <- pheno[[col]]
    
    # Skip if mostly NA
    if (sum(is.na(col_data)) > nrow(pheno) * 0.5) next
    
    # Check if categorical with reasonable number of levels
    unique_vals <- unique(col_data[!is.na(col_data)])
    n_unique <- length(unique_vals)
    
    if (n_unique >= 2 && n_unique <= min(20, nrow(pheno) / 2)) {
      grouping_vars[[col]] <- list(
        n_levels = n_unique,
        levels = as.character(unique_vals),
        counts = table(col_data)
      )
    }
  }
  
  return(grouping_vars)
}

#' Get SRA accessions linked to GEO samples
#' 
#' Extracts SRA (Sequence Read Archive) accession numbers from GEO metadata.
#' Useful for downloading raw sequencing data associated with GEO entries.
#' 
#' @param geo_object GEO object (ExpressionSet or list) or accession string
#' @param return_links Logical, return full SRA links instead of just accessions (default: FALSE)
#' @return Data frame with sample IDs and corresponding SRA accessions
#' @export
#' @examples
#' \dontrun{
#' gse <- GEOquery::getGEO("GSE12345")
#' sra_info <- get_geo_sra_links(gse)
#' 
#' # Or directly from accession
#' sra_info <- get_geo_sra_links("GSE12345")
#' }
get_geo_sra_links <- function(geo_object, return_links = FALSE) {
  
  # If string provided, download the object first
  if (is.character(geo_object)) {
    message("Downloading GEO data for ", geo_object, "...")
    geo_object <- GEOquery::getGEO(geo_object)
  }
  
  # Extract phenotype data
  if (inherits(geo_object, "list")) {
    pheno <- Biobase::pData(geo_object[[1]])
  } else if (inherits(geo_object, "ExpressionSet")) {
    pheno <- Biobase::pData(geo_object)
  } else {
    stop("Input must be a GEO object or GEO accession string")
  }
  
  # Look for SRA accessions in various columns
  sra_columns <- grep("sra|run|experiment|submission", 
                     colnames(pheno), value = TRUE, ignore.case = TRUE)
  
  result <- data.frame(
    sample_id = rownames(pheno),
    stringsAsFactors = FALSE
  )
  
  # Extract SRA accessions
  sra_accessions <- character(nrow(pheno))
  
  for (col in sra_columns) {
    col_data <- as.character(pheno[[col]])
    
    # Extract SRX, SRR, SRS, or SRP accessions
    for (i in seq_along(col_data)) {
      if (is.na(sra_accessions[i]) || sra_accessions[i] == "") {
        matches <- regmatches(col_data[i], 
                            gregexpr("SR[XRSP]\\d+", col_data[i]))[[1]]
        if (length(matches) > 0) {
          sra_accessions[i] <- paste(matches, collapse = ";")
        }
      }
    }
  }
  
  result$sra_accession <- sra_accessions
  
  # Add links if requested
  if (return_links) {
    result$sra_link <- ifelse(
      result$sra_accession != "",
      paste0("https://www.ncbi.nlm.nih.gov/sra/", result$sra_accession),
      NA
    )
  }
  
  # Filter to only samples with SRA data
  result <- result[result$sra_accession != "", ]
  
  if (nrow(result) == 0) {
    message("No SRA accessions found in GEO metadata")
  } else {
    message("Found SRA accessions for ", nrow(result), " sample(s)")
  }
  
  return(result)
}

#' Search GEO for datasets matching criteria
#' 
#' Search GEO database using various criteria including organism, tissue,
#' experiment type, and keywords.
#' 
#' @param term Search term or keywords
#' @param organism Organism name (e.g., "Homo sapiens", "Mus musculus")
#' @param entry_type Type of GEO entry: "gse", "gds", "gpl", "gsm" (default: "gse")
#' @param max_results Maximum number of results to return (default: 100)
#' @return Data frame with search results
#' @export
#' @examples
#' \dontrun{
#' # Search for human RNA-seq datasets
#' results <- search_geo(term = "RNA-seq", organism = "Homo sapiens")
#' 
#' # Search for cancer datasets
#' results <- search_geo(term = "cancer", max_results = 50)
#' }
search_geo <- function(term, 
                      organism = NULL,
                      entry_type = c("gse", "gds", "gpl", "gsm"),
                      max_results = 100) {
  
  entry_type <- match.arg(entry_type)
  
  # Build search query
  query <- term
  if (!is.null(organism)) {
    query <- paste0(query, " AND ", organism, "[Organism]")
  }
  query <- paste0(query, " AND ", entry_type, "[Entry Type]")
  
  message("Searching GEO for: ", query)
  message("This requires rentrez package...")
  
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    stop("Package 'rentrez' is required for GEO search. Install it with: install.packages('rentrez')")
  }
  
  tryCatch({
    # Search using NCBI E-utilities
    search_results <- rentrez::entrez_search(
      db = "gds",
      term = query,
      retmax = max_results
    )
    
    if (length(search_results$ids) == 0) {
      message("No results found")
      return(data.frame())
    }
    
    message("Found ", length(search_results$ids), " result(s)")
    message("Fetching details...")
    
    # Fetch summaries
    summaries <- rentrez::entrez_summary(db = "gds", id = search_results$ids)
    
    # Parse results
    results <- lapply(summaries, function(x) {
      data.frame(
        accession = x$accession,
        title = x$title,
        summary = x$summary,
        entry_type = x$entrytype,
        organism = paste(x$taxon, collapse = "; "),
        platform = x$gpl,
        n_samples = ifelse(!is.null(x$n_samples), x$n_samples, NA),
        pubmed_id = ifelse(!is.null(x$pubmedids), paste(x$pubmedids, collapse = "; "), NA),
        stringsAsFactors = FALSE
      )
    })
    
    results_df <- do.call(rbind, results)
    return(results_df)
    
  }, error = function(e) {
    stop("GEO search failed: ", e$message)
  })
}

#' Convert GEO expression data to different formats
#' 
#' Convert GEO ExpressionSet to various formats including data.frame,
#' SummarizedExperiment, or files (CSV, TSV, RDS).
#' 
#' @param geo_object GEO object (ExpressionSet or list)
#' @param format Output format: "data.frame", "matrix", "SummarizedExperiment", 
#'   "csv", "tsv", "rds" (default: "data.frame")
#' @param output_file File path for file formats (csv, tsv, rds)
#' @param include_pheno Logical, include phenotype data (default: TRUE)
#' @param log_transform Logical, log2 transform expression values (default: FALSE)
#' @return Converted data object or file path
#' @export
#' @examples
#' \dontrun{
#' gse <- GEOquery::getGEO("GSE12345")
#' 
#' # Convert to data frame
#' df <- convert_geo_format(gse, format = "data.frame")
#' 
#' # Save as CSV
#' convert_geo_format(gse, format = "csv", output_file = "expression.csv")
#' 
#' # Convert to SummarizedExperiment
#' se <- convert_geo_format(gse, format = "SummarizedExperiment")
#' }
convert_geo_format <- function(geo_object,
                              format = c("data.frame", "matrix", "SummarizedExperiment", 
                                        "csv", "tsv", "rds"),
                              output_file = NULL,
                              include_pheno = TRUE,
                              log_transform = FALSE) {
  
  format <- match.arg(format)
  
  # Extract ExpressionSet
  if (inherits(geo_object, "list")) {
    if (length(geo_object) > 1) {
      warning("Multiple platforms detected. Using first platform: ", 
             Biobase::annotation(geo_object[[1]]))
    }
    eset <- geo_object[[1]]
  } else if (inherits(geo_object, "ExpressionSet")) {
    eset <- geo_object
  } else {
    stop("Input must be ExpressionSet or list of ExpressionSets")
  }
  
  # Get expression matrix
  expr <- Biobase::exprs(eset)
  
  # Log transform if requested
  if (log_transform) {
    message("Applying log2 transformation...")
    expr <- log2(expr + 1)
  }
  
  # Get phenotype data
  pheno <- if (include_pheno) Biobase::pData(eset) else NULL
  
  # Convert based on format
  result <- switch(format,
    "matrix" = {
      expr
    },
    "data.frame" = {
      df <- as.data.frame(expr)
      if (include_pheno) {
        attr(df, "phenotype_data") <- pheno
      }
      df
    },
    "SummarizedExperiment" = {
      if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("Package 'SummarizedExperiment' required. Install from Bioconductor.")
      }
      
      fdata <- Biobase::fData(eset)
      row_data <- if (nrow(fdata) > 0) fdata else NULL
      
      SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = expr),
        colData = pheno,
        rowData = row_data
      )
    },
    "csv" = {
      if (is.null(output_file)) {
        output_file <- paste0(Biobase::annotation(eset), "_expression.csv")
      }
      write.csv(expr, output_file, row.names = TRUE)
      if (include_pheno) {
        pheno_file <- gsub("\\.csv$", "_phenotype.csv", output_file)
        write.csv(pheno, pheno_file, row.names = TRUE)
        message("Phenotype data saved to: ", pheno_file)
      }
      message("Expression data saved to: ", output_file)
      output_file
    },
    "tsv" = {
      if (is.null(output_file)) {
        output_file <- paste0(Biobase::annotation(eset), "_expression.tsv")
      }
      write.table(expr, output_file, sep = "\t", row.names = TRUE, quote = FALSE)
      if (include_pheno) {
        pheno_file <- gsub("\\.tsv$", "_phenotype.tsv", output_file)
        write.table(pheno, pheno_file, sep = "\t", row.names = TRUE, quote = FALSE)
        message("Phenotype data saved to: ", pheno_file)
      }
      message("Expression data saved to: ", output_file)
      output_file
    },
    "rds" = {
      if (is.null(output_file)) {
        output_file <- paste0(Biobase::annotation(eset), "_data.rds")
      }
      data_list <- list(
        expression = expr,
        phenotype = pheno,
        platform = Biobase::annotation(eset)
      )
      saveRDS(data_list, output_file)
      message("Data saved to: ", output_file)
      output_file
    }
  )
  
  return(result)
}

#' Batch download multiple GEO datasets
#' 
#' Download multiple GEO datasets in batch with progress tracking and error handling.
#' 
#' @param accessions Character vector of GEO accessions
#' @param output_dir Base output directory (default: "geo_batch")
#' @param data_type Type of data to download (see download_geo)
#' @param parallel Logical, download in parallel (default: FALSE)
#' @param n_cores Number of cores for parallel processing (default: 2)
#' @param continue_on_error Logical, continue if one download fails (default: TRUE)
#' @param ... Additional arguments passed to download_geo
#' @return List of download results for each accession
#' @export
#' @examples
#' \dontrun{
#' # Download multiple datasets
#' accessions <- c("GSE12345", "GSE67890", "GSE11111")
#' results <- batch_download_geo(accessions, data_type = "matrix")
#' 
#' # Check which succeeded
#' succeeded <- sapply(results, function(x) x$success)
#' }
batch_download_geo <- function(accessions,
                              output_dir = "geo_batch",
                              data_type = "matrix",
                              parallel = FALSE,
                              n_cores = 2,
                              continue_on_error = TRUE,
                              ...) {
  
  check_and_create_dir(output_dir)
  
  n_accessions <- length(accessions)
  message("Starting batch download of ", n_accessions, " GEO dataset(s)")
  message(paste0(strrep("=", 70)))
  
  results <- list()
  
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' not available. Using sequential download.")
      parallel <- FALSE
    } else {
      message("Using parallel processing with ", n_cores, " cores")
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl))
      
      results <- parallel::parLapply(cl, accessions, function(acc) {
        tryCatch({
          download_geo(acc, output_dir = file.path(output_dir, acc), 
                      data_type = data_type, verbose = FALSE, ...)
        }, error = function(e) {
          list(accession = acc, success = FALSE, error = e$message)
        })
      })
      names(results) <- accessions
    }
  }
  
  if (!parallel) {
    for (i in seq_along(accessions)) {
      acc <- accessions[i]
      message("\n[", i, "/", n_accessions, "] Processing: ", acc)
      
      result <- tryCatch({
        download_geo(acc, output_dir = file.path(output_dir, acc),
                    data_type = data_type, ...)
      }, error = function(e) {
        message("ERROR: ", e$message)
        if (!continue_on_error) {
          stop("Batch download stopped due to error")
        }
        list(accession = acc, success = FALSE, error = e$message)
      })
      
      results[[acc]] <- result
    }
  }
  
  # Summary
  n_success <- sum(sapply(results, function(x) isTRUE(x$success)))
  message("\n", paste0(strrep("=", 70)))
  message("BATCH DOWNLOAD COMPLETE")
  message("Successful: ", n_success, "/", n_accessions)
  message("Failed: ", n_accessions - n_success, "/", n_accessions)
  message(paste0(strrep("=", 70)))
  
  return(results)
}
