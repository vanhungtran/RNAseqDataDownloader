#' Download data from TCGA
#' 
#' Advanced function to download various types of TCGA data including RNA-seq,
#' miRNA-seq, methylation, copy number, mutations, and clinical data. Supports
#' multiple workflows, sample type filtering, and comprehensive error handling.
#' 
#' @param project TCGA project ID (e.g., "TCGA-BRCA", "TCGA-LUAD"). Use list_tcga_projects() 
#'   to see available projects
#' @param output_dir Directory to save downloaded data (default: "tcga_data")
#' @param data_category Data category to download. Options:
#'   \itemize{
#'     \item "Transcriptome Profiling" - RNA-seq data (default)
#'     \item "Simple Nucleotide Variation" - Mutation data
#'     \item "Copy Number Variation" - CNV data
#'     \item "DNA Methylation" - Methylation arrays
#'     \item "Clinical" - Clinical/phenotype data
#'     \item "Biospecimen" - Sample metadata
#'   }
#' @param data_type Specific data type within category. Common options:
#'   \itemize{
#'     \item "Gene Expression Quantification" (default for RNA-seq)
#'     \item "Isoform Expression Quantification"
#'     \item "miRNA Expression Quantification"
#'     \item "Masked Somatic Mutation"
#'     \item "Gene Level Copy Number"
#'     \item "Methylation Beta Value"
#'   }
#' @param workflow_type Analysis workflow. Options:
#'   \itemize{
#'     \item "STAR - Counts" - STAR aligned counts (recommended)
#'     \item "HTSeq - Counts" - HTSeq counts
#'     \item "HTSeq - FPKM" - FPKM values
#'     \item "HTSeq - FPKM-UQ" - Upper quartile normalized FPKM
#'   }
#' @param sample_type Sample type filter. Options: "Primary Tumor", "Solid Tissue Normal",
#'   "Metastatic", "Blood Derived Normal", etc. Use NULL for all types
#' @param experimental_strategy Experimental strategy: "RNA-Seq", "miRNA-Seq", 
#'   "Genotyping Array", "Methylation Array", etc.
#' @param barcode Specific sample barcodes to download (character vector)
#' @param legacy Use legacy TCGA data archive (default: FALSE for harmonized data)
#' @param download_clinical Also download clinical data (default: TRUE)
#' @param download_biospecimen Also download biospecimen data (default: FALSE)
#' @param prepare_data Prepare/process downloaded data (default: TRUE)
#' @param save_format Save format: "rds", "csv", "tsv", or "all" (default: "rds")
#' @param force_download Re-download even if files exist (default: FALSE)
#' @param max_retries Maximum download retry attempts (default: 3)
#' @param verbose Print detailed progress messages (default: TRUE)
#' @param ... Additional parameters passed to GDCquery
#' @return List containing:
#'   \itemize{
#'     \item data - Prepared SummarizedExperiment object (if prepare_data = TRUE)
#'     \item query - GDC query object
#'     \item files - Downloaded file paths
#'     \item clinical - Clinical data (if requested)
#'     \item metadata - Summary metadata
#'   }
#' @export
#' @examples
#' \dontrun{
#' # Download BRCA RNA-seq data (STAR counts)
#' result <- download_tcga("TCGA-BRCA", workflow_type = "STAR - Counts")
#' 
#' # Download only tumor samples
#' result <- download_tcga("TCGA-LUAD", sample_type = "Primary Tumor")
#' 
#' # Download mutation data
#' result <- download_tcga("TCGA-BRCA", 
#'                        data_category = "Simple Nucleotide Variation",
#'                        data_type = "Masked Somatic Mutation")
#' 
#' # Download methylation data
#' result <- download_tcga("TCGA-BRCA",
#'                        data_category = "DNA Methylation",
#'                        experimental_strategy = "Methylation Array")
#' }
download_tcga <- function(project, 
                         output_dir = "tcga_data",
                         data_category = "Transcriptome Profiling",
                         data_type = "Gene Expression Quantification",
                         workflow_type = "STAR - Counts",
                         sample_type = NULL,
                         experimental_strategy = NULL,
                         barcode = NULL,
                         legacy = FALSE,
                         download_clinical = TRUE,
                         download_biospecimen = FALSE,
                         prepare_data = TRUE,
                         save_format = c("rds", "csv", "tsv", "all"),
                         force_download = FALSE,
                         max_retries = 3,
                         verbose = TRUE,
                         ...) {
  
  save_format <- match.arg(save_format)
  
  # Validate project ID
  project <- toupper(trimws(project))
  if (!grepl("^TCGA-[A-Z]+$", project)) {
    stop("Invalid TCGA project ID format. Should be 'TCGA-XXX' (e.g., 'TCGA-BRCA')")
  }
  
  # Create project-specific directory
  check_and_create_dir(output_dir)
  project_dir <- file.path(output_dir, project)
  check_and_create_dir(project_dir)
  
  if (verbose) {
    message(paste0(strrep("=", 70)))
    message("TCGA Data Download")
    message("Project: ", project)
    message("Category: ", data_category)
    message("Type: ", data_type)
    if (!is.null(workflow_type)) message("Workflow: ", workflow_type)
    if (!is.null(sample_type)) message("Sample Type: ", paste(sample_type, collapse = ", "))
    message("Output: ", project_dir)
    message(paste0(strrep("=", 70)))
  }
  
  # Initialize result
  result <- list(
    project = project,
    data_category = data_category,
    data_type = data_type,
    query = NULL,
    data = NULL,
    clinical = NULL,
    biospecimen = NULL,
    files = character(0),
    metadata = list(),
    success = FALSE
  )
  
  tryCatch({
    
    # Build query
    if (verbose) message("\n[1/", ifelse(prepare_data, "4", "3"), "] Building GDC query...")
    
    query_args <- list(
      project = project,
      data.category = data_category,
      data.type = data_type,
      legacy = legacy
    )
    
    if (!is.null(workflow_type)) query_args$workflow.type <- workflow_type
    if (!is.null(sample_type)) query_args$sample.type <- sample_type
    if (!is.null(experimental_strategy)) query_args$experimental.strategy <- experimental_strategy
    if (!is.null(barcode)) query_args$barcode <- barcode
    
    # Add additional arguments
    query_args <- c(query_args, list(...))
    
    query <- do.call(TCGAbiolinks::GDCquery, query_args)
    result$query <- query
    
    # Get query results summary
    query_results <- TCGAbiolinks::getResults(query)
    
    if (verbose) {
      message("Query returned ", nrow(query_results), " file(s)")
      if (nrow(query_results) > 0) {
        message("Total size: ", format(sum(query_results$file_size) / 1e9, digits = 3), " GB")
        
        # Show sample type breakdown
        if ("sample_type" %in% colnames(query_results)) {
          sample_counts <- table(query_results$sample_type)
          message("Sample types:")
          for (st in names(sample_counts)) {
            message("  - ", st, ": ", sample_counts[st])
          }
        }
      }
    }
    
    if (nrow(query_results) == 0) {
      warning("No data found matching the query criteria")
      return(result)
    }
    
    # Download data with retry logic
    if (verbose) message("\n[2/", ifelse(prepare_data, "4", "3"), "] Downloading data from GDC...")
    
    attempt <- 0
    download_success <- FALSE
    
    while (attempt < max_retries && !download_success) {
      attempt <- attempt + 1
      if (verbose && attempt > 1) message("Retry attempt ", attempt, "/", max_retries)
      
      tryCatch({
        TCGAbiolinks::GDCdownload(
          query, 
          directory = project_dir,
          method = "api",
          files.per.chunk = 10
        )
        download_success <- TRUE
        if (verbose) message("Download completed successfully")
      }, error = function(e) {
        if (attempt >= max_retries) {
          stop("Failed to download after ", max_retries, " attempts: ", e$message)
        }
        if (verbose) message("Download failed: ", e$message)
        Sys.sleep(2 ^ attempt)
      })
    }
    
    # List downloaded files
    gdc_dir <- file.path(project_dir, "GDCdata", project)
    if (dir.exists(gdc_dir)) {
      result$files <- list.files(gdc_dir, recursive = TRUE, full.names = TRUE)
    }
    
    # Prepare data
    if (prepare_data) {
      if (verbose) message("\n[3/4] Preparing and processing data...")
      
      tryCatch({
        data <- TCGAbiolinks::GDCprepare(
          query, 
          directory = project_dir,
          summarizedExperiment = TRUE
        )
        result$data <- data
        
        # Extract metadata
        result$metadata <- .extract_tcga_metadata(data, verbose)
        
        # Save prepared data
        .save_tcga_data(data, project_dir, project, save_format, verbose)
        
        if (verbose) message("Data preparation completed")
        
      }, error = function(e) {
        warning("Failed to prepare data: ", e$message)
      })
    }
    
    # Download clinical data
    if (download_clinical) {
      if (verbose) message("\n[", ifelse(prepare_data, "4", "3"), "/", 
                          ifelse(prepare_data, "4", "3"), "] Downloading clinical data...")
      
      tryCatch({
        clinical <- TCGAbiolinks::GDCquery_clinic(project, type = "clinical")
        result$clinical <- clinical
        
        # Save clinical data
        clinical_file <- file.path(project_dir, paste0(project, "_clinical.tsv"))
        write.table(clinical, clinical_file, sep = "\t", row.names = FALSE, quote = FALSE)
        
        if (verbose) message("Clinical data: ", nrow(clinical), " patients")
        
      }, error = function(e) {
        warning("Failed to download clinical data: ", e$message)
      })
    }
    
    # Download biospecimen data
    if (download_biospecimen) {
      if (verbose) message("\nDownloading biospecimen data...")
      
      tryCatch({
        biospecimen <- TCGAbiolinks::GDCquery_clinic(project, type = "biospecimen")
        result$biospecimen <- biospecimen
        
        # Save biospecimen data
        biospecimen_file <- file.path(project_dir, paste0(project, "_biospecimen.tsv"))
        write.table(biospecimen, biospecimen_file, sep = "\t", row.names = FALSE, quote = FALSE)
        
        if (verbose) message("Biospecimen data: ", nrow(biospecimen), " samples")
        
      }, error = function(e) {
        warning("Failed to download biospecimen data: ", e$message)
      })
    }
    
    result$success <- TRUE
    
    if (verbose) {
      message("\n", paste0(strrep("=", 70)))
      message("SUCCESS: TCGA data download completed")
      message("Project: ", project)
      message("Files saved to: ", project_dir)
      message(paste0(strrep("=", 70)))
    }
    
    return(result)
    
  }, error = function(e) {
    warning("Failed to download TCGA data: ", e$message)
    result$error <- e$message
    return(result)
  })
}

#' List available TCGA projects
#' 
#' Get a list of all available TCGA projects with descriptions.
#' 
#' @param include_count Include number of cases per project (default: TRUE)
#' @param legacy Use legacy archive (default: FALSE)
#' @return Data frame with project information
#' @export
#' @examples
#' \dontrun{
#' # List all TCGA projects
#' projects <- list_tcga_projects()
#' 
#' # View cancer types
#' View(projects)
#' }
list_tcga_projects <- function(include_count = TRUE, legacy = FALSE) {
  
  message("Fetching TCGA project list from GDC...")
  
  tryCatch({
    projects <- TCGAbiolinks::getGDCprojects()
    
    # Filter for TCGA projects
    tcga_projects <- projects[grepl("^TCGA-", projects$project_id), ]
    
    if (include_count) {
      # Get case counts
      tcga_projects$n_cases <- sapply(tcga_projects$project_id, function(proj) {
        tryCatch({
          query <- TCGAbiolinks::GDCquery(
            project = proj,
            data.category = "Transcriptome Profiling",
            data.type = "Gene Expression Quantification",
            legacy = legacy
          )
          results <- TCGAbiolinks::getResults(query)
          length(unique(results$cases))
        }, error = function(e) NA)
      })
    }
    
    # Sort by project ID
    tcga_projects <- tcga_projects[order(tcga_projects$project_id), ]
    
    message("Found ", nrow(tcga_projects), " TCGA project(s)")
    
    return(tcga_projects)
    
  }, error = function(e) {
    stop("Failed to fetch TCGA projects: ", e$message)
  })
}

#' Search TCGA for samples matching criteria
#' 
#' Search TCGA database for samples based on various clinical and molecular criteria.
#' 
#' @param project TCGA project ID
#' @param gender Patient gender: "male", "female", or NULL for all
#' @param race Patient race filter
#' @param vital_status "Alive" or "Dead"
#' @param age_range Age range as c(min, max)
#' @param tumor_stage Tumor stage (e.g., "Stage I", "Stage II", etc.)
#' @param sample_type Sample type filter
#' @return Data frame with matching samples and clinical information
#' @export
#' @examples
#' \dontrun{
#' # Search for female BRCA patients
#' samples <- search_tcga_samples("TCGA-BRCA", gender = "female")
#' 
#' # Search for Stage I patients
#' samples <- search_tcga_samples("TCGA-LUAD", tumor_stage = "Stage I")
#' }
search_tcga_samples <- function(project,
                               gender = NULL,
                               race = NULL,
                               vital_status = NULL,
                               age_range = NULL,
                               tumor_stage = NULL,
                               sample_type = NULL) {
  
  message("Searching TCGA samples for ", project)
  
  tryCatch({
    # Get clinical data
    clinical <- TCGAbiolinks::GDCquery_clinic(project, type = "clinical")
    
    # Apply filters
    if (!is.null(gender)) {
      clinical <- clinical[tolower(clinical$gender) == tolower(gender), ]
    }
    
    if (!is.null(race)) {
      clinical <- clinical[grepl(race, clinical$race, ignore.case = TRUE), ]
    }
    
    if (!is.null(vital_status)) {
      clinical <- clinical[tolower(clinical$vital_status) == tolower(vital_status), ]
    }
    
    if (!is.null(age_range) && length(age_range) == 2) {
      age_col <- if("age_at_diagnosis" %in% colnames(clinical)) "age_at_diagnosis" else "age_at_index"
      if (age_col %in% colnames(clinical)) {
        clinical <- clinical[clinical[[age_col]] >= age_range[1] & 
                           clinical[[age_col]] <= age_range[2], ]
      }
    }
    
    if (!is.null(tumor_stage)) {
      stage_cols <- grep("stage|tumor", colnames(clinical), ignore.case = TRUE, value = TRUE)
      if (length(stage_cols) > 0) {
        matched <- apply(clinical[, stage_cols, drop = FALSE], 1, function(x) {
          any(grepl(tumor_stage, x, ignore.case = TRUE))
        })
        clinical <- clinical[matched, ]
      }
    }
    
    message("Found ", nrow(clinical), " matching sample(s)")
    
    return(clinical)
    
  }, error = function(e) {
    stop("Search failed: ", e$message)
  })
}

#' Extract comprehensive metadata from TCGA SummarizedExperiment
#' @keywords internal
.extract_tcga_metadata <- function(se_object, verbose = TRUE) {
  
  metadata <- list()
  
  # Basic information
  metadata$n_samples <- ncol(se_object)
  metadata$n_features <- nrow(se_object)
  metadata$sample_ids <- colnames(se_object)
  metadata$feature_ids <- head(rownames(se_object), 100)
  
  # Column data (clinical/sample info)
  if (!is.null(SummarizedExperiment::colData(se_object))) {
    coldata <- as.data.frame(SummarizedExperiment::colData(se_object))
    metadata$phenotype_variables <- colnames(coldata)
    metadata$n_phenotype_vars <- ncol(coldata)
    
    # Sample type distribution
    if ("sample_type" %in% colnames(coldata)) {
      metadata$sample_types <- table(coldata$sample_type)
    }
    
    # Tumor stage distribution
    stage_cols <- grep("stage", colnames(coldata), ignore.case = TRUE, value = TRUE)
    if (length(stage_cols) > 0) {
      metadata$tumor_stages <- table(coldata[[stage_cols[1]]])
    }
  }
  
  # Row data (feature annotation)
  if (!is.null(SummarizedExperiment::rowData(se_object))) {
    rowdata <- as.data.frame(SummarizedExperiment::rowData(se_object))
    metadata$feature_annotation_vars <- colnames(rowdata)
  }
  
  # Assay information
  metadata$assay_names <- SummarizedExperiment::assayNames(se_object)
  
  # Expression statistics
  if (length(metadata$assay_names) > 0) {
    expr <- SummarizedExperiment::assay(se_object, 1)
    metadata$expression_stats <- list(
      min = min(expr, na.rm = TRUE),
      max = max(expr, na.rm = TRUE),
      mean = mean(expr, na.rm = TRUE),
      median = median(expr, na.rm = TRUE),
      na_count = sum(is.na(expr)),
      zero_count = sum(expr == 0, na.rm = TRUE)
    )
  }
  
  if (verbose) {
    message("Samples: ", metadata$n_samples)
    message("Features: ", metadata$n_features)
    if (!is.null(metadata$sample_types)) {
      message("Sample types: ", paste(names(metadata$sample_types), collapse = ", "))
    }
  }
  
  return(metadata)
}

#' Save TCGA data in various formats
#' @keywords internal
.save_tcga_data <- function(se_object, output_dir, project, format, verbose) {
  
  formats_to_save <- if (format == "all") c("rds", "csv", "tsv") else format
  
  for (fmt in formats_to_save) {
    tryCatch({
      if (fmt == "rds") {
        # Save as RDS
        rds_file <- file.path(output_dir, paste0(project, "_data.rds"))
        saveRDS(se_object, rds_file)
        if (verbose) message("Saved RDS: ", basename(rds_file))
        
      } else if (fmt %in% c("csv", "tsv")) {
        # Extract expression matrix
        expr <- SummarizedExperiment::assay(se_object, 1)
        coldata <- as.data.frame(SummarizedExperiment::colData(se_object))
        
        sep <- if (fmt == "csv") "," else "\t"
        ext <- if (fmt == "csv") ".csv" else ".tsv"
        
        # Save expression data
        expr_file <- file.path(output_dir, paste0(project, "_expression", ext))
        write.table(expr, expr_file, sep = sep, row.names = TRUE, 
                   col.names = NA, quote = FALSE)
        
        # Save sample data
        sample_file <- file.path(output_dir, paste0(project, "_samples", ext))
        write.table(coldata, sample_file, sep = sep, row.names = TRUE, 
                   col.names = NA, quote = FALSE)
        
        if (verbose) message("Saved ", toupper(fmt), ": ", basename(expr_file), 
                           " and ", basename(sample_file))
      }
    }, error = function(e) {
      warning("Failed to save ", fmt, " format: ", e$message)
    })
  }
}

#' Batch download multiple TCGA projects
#' 
#' Download data from multiple TCGA projects in batch with progress tracking.
#' 
#' @param projects Character vector of TCGA project IDs
#' @param output_dir Base output directory (default: "tcga_batch")
#' @param parallel Process projects in parallel (default: FALSE)
#' @param n_cores Number of cores for parallel processing (default: 2)
#' @param continue_on_error Continue if one project fails (default: TRUE)
#' @param ... Additional arguments passed to download_tcga
#' @return List of download results for each project
#' @export
#' @examples
#' \dontrun{
#' # Download multiple cancer types
#' projects <- c("TCGA-BRCA", "TCGA-LUAD", "TCGA-PRAD")
#' results <- batch_download_tcga(projects)
#' 
#' # Check success
#' sapply(results, function(x) x$success)
#' }
batch_download_tcga <- function(projects,
                               output_dir = "tcga_batch",
                               parallel = FALSE,
                               n_cores = 2,
                               continue_on_error = TRUE,
                               ...) {
  
  check_and_create_dir(output_dir)
  
  n_projects <- length(projects)
  message("Starting batch download of ", n_projects, " TCGA project(s)")
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
      
      results <- parallel::parLapply(cl, projects, function(proj) {
        tryCatch({
          download_tcga(proj, output_dir = output_dir, verbose = FALSE, ...)
        }, error = function(e) {
          list(project = proj, success = FALSE, error = e$message)
        })
      })
      names(results) <- projects
    }
  }
  
  if (!parallel) {
    for (i in seq_along(projects)) {
      proj <- projects[i]
      message("\n[", i, "/", n_projects, "] Processing: ", proj)
      
      result <- tryCatch({
        download_tcga(proj, output_dir = output_dir, ...)
      }, error = function(e) {
        message("ERROR: ", e$message)
        if (!continue_on_error) {
          stop("Batch download stopped due to error")
        }
        list(project = proj, success = FALSE, error = e$message)
      })
      
      results[[proj]] <- result
    }
  }
  
  # Summary
  n_success <- sum(sapply(results, function(x) isTRUE(x$success)))
  message("\n", paste0(strrep("=", 70)))
  message("BATCH DOWNLOAD COMPLETE")
  message("Successful: ", n_success, "/", n_projects)
  message("Failed: ", n_projects - n_success, "/", n_projects)
  message(paste0(strrep("=", 70)))
  
  return(results)
}
