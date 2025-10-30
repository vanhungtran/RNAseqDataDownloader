#' Download data from Single Cell Expression Atlas
#' @param accession Experiment accession number (e.g., "E-MTAB-5061")
#' @param output_dir Directory to save downloaded data
#' @param file_type Type of file to download: "loom", "rds", or "processed"
#' @return Path to downloaded files
#' @export
download_sc_atlas <- function(accession, output_dir = "sc_data", file_type = "loom") {
  check_and_create_dir(output_dir)

  tryCatch({
    # Check file type
    if (!file_type %in% c("loom", "rds", "processed")) {
      stop("file_type must be one of: loom, rds, processed")
    }

    # Construct URL
    base_url <- "https://www.ebi.ac.uk/gxa/sc/experiment/"
    url <- paste0(base_url, accession, "/download/", file_type)

    # Download file
    dest_file <- file.path(output_dir, paste0(accession, ".", file_type))
    download.file(url, destfile = dest_file)

    message("Single Cell Atlas data downloaded successfully")
    return(dest_file)
  }, error = function(e) {
    stop("Failed to download Single Cell Atlas data: ", e$message)
  })
}

#' Download data from Human Cell Atlas (CZ CELLxGENE)
#' 
#' Downloads single-cell data from the CZ CELLxGENE Discover portal, which hosts
#' Human Cell Atlas and other single-cell datasets. This function uses the modern
#' cellxgenedp package (replacement for obsolete cellexalvrR).
#' 
#' @param dataset_id CELLxGENE dataset ID (e.g., "b9fc3d70-5a72-4479-a046-c2cc1ab19efc")
#'   Use `cellxgenedp::datasets()` to browse available datasets
#' @param output_dir Directory to save downloaded data (default: "hca_data")
#' @param file_format File format: "h5ad" (default, AnnData) or "rds" (Seurat/SingleCellExperiment)
#' @param organism Filter by organism: "Homo sapiens", "Mus musculus", or NULL for all
#' @param tissue Filter by tissue type, or NULL for all
#' @param disease Filter by disease, or NULL for all
#' @param verbose Print progress messages (default: TRUE)
#' @return List containing:
#'   \itemize{
#'     \item data_path - Path to downloaded file
#'     \item dataset_info - Dataset metadata
#'     \item file_path - Full path to downloaded file
#'   }
#' @export
#' @examples
#' \dontrun{
#' # Browse available datasets
#' library(cellxgenedp)
#' datasets <- datasets()
#' head(datasets)
#' 
#' # Download a specific dataset
#' result <- download_hca(
#'   dataset_id = "b9fc3d70-5a72-4479-a046-c2cc1ab19efc",
#'   output_dir = "hca_data"
#' )
#' 
#' # Load the downloaded data
#' library(Seurat)
#' adata <- Read10X_h5(result$file_path)
#' seurat_obj <- CreateSeuratObject(counts = adata)
#' }
download_hca <- function(dataset_id = NULL, 
                        output_dir = "hca_data", 
                        file_format = c("h5ad", "rds"),
                        organism = NULL,
                        tissue = NULL,
                        disease = NULL,
                        verbose = TRUE) {
  
  # Check if cellxgenedp is available
  if (!requireNamespace("cellxgenedp", quietly = TRUE)) {
    stop("Package 'cellxgenedp' is required for HCA downloads.\n",
         "Install it with: BiocManager::install('cellxgenedp')")
  }
  
  file_format <- match.arg(file_format)
  check_and_create_dir(output_dir)
  
  tryCatch({
    # Get available datasets
    if (verbose) {
      message("Fetching CELLxGENE datasets catalog...")
    }
    
    all_datasets <- cellxgenedp::datasets()
    
    # Apply filters if specified
    filtered_datasets <- all_datasets
    
    if (!is.null(organism)) {
      filtered_datasets <- filtered_datasets[
        grepl(organism, filtered_datasets$organism, ignore.case = TRUE), 
      ]
    }
    
    if (!is.null(tissue)) {
      filtered_datasets <- filtered_datasets[
        grepl(tissue, filtered_datasets$tissue, ignore.case = TRUE), 
      ]
    }
    
    if (!is.null(disease)) {
      filtered_datasets <- filtered_datasets[
        grepl(disease, filtered_datasets$disease, ignore.case = TRUE), 
      ]
    }
    
    # If no dataset_id provided, show available datasets
    if (is.null(dataset_id)) {
      if (verbose) {
        message("No dataset_id provided. Showing available datasets:")
        message("Total datasets found: ", nrow(filtered_datasets))
        print(head(filtered_datasets[, c("dataset_id", "dataset_title", 
                                        "organism", "tissue")], 10))
      }
      return(list(
        available_datasets = filtered_datasets,
        message = "Please specify a dataset_id to download"
      ))
    }
    
    # Validate dataset_id
    if (!dataset_id %in% all_datasets$dataset_id) {
      stop("Dataset ID '", dataset_id, "' not found in CELLxGENE catalog")
    }
    
    # Get dataset info
    dataset_info <- all_datasets[all_datasets$dataset_id == dataset_id, ]
    
    if (verbose) {
      message("Dataset: ", dataset_info$dataset_title)
      message("Organism: ", dataset_info$organism)
      message("Tissue: ", dataset_info$tissue)
      message("Cell count: ", dataset_info$cell_count)
    }
    
    # Download the dataset
    if (verbose) message("Downloading dataset...")
    
    # Get files for this dataset
    files_db <- cellxgenedp::files()
    dataset_files <- files_db[files_db$dataset_id == dataset_id, ]
    
    if (nrow(dataset_files) == 0) {
      stop("No files found for dataset: ", dataset_id)
    }
    
    # Filter by format
    if (file_format == "h5ad") {
      target_file <- dataset_files[dataset_files$filetype == "H5AD", ]
    } else if (file_format == "rds") {
      target_file <- dataset_files[dataset_files$filetype == "RDS", ]
    }
    
    if (nrow(target_file) == 0) {
      warning("Format '", file_format, "' not available. Using first available file.")
      target_file <- dataset_files[1, ]
    } else {
      target_file <- target_file[1, ]  # Use first matching file
    }
    
    # Download the file
    local_file <- cellxgenedp::files_download(
      target_file$file_id,
      destdir = output_dir
    )
    
    if (verbose) {
      message("Download complete!")
      message("File saved to: ", local_file)
    }
    
    result <- list(
      data_path = output_dir,
      file_path = local_file,
      dataset_info = as.list(dataset_info),
      file_format = file_format,
      success = TRUE
    )
    
    return(result)
    
  }, error = function(e) {
    message("Error downloading HCA data: ", e$message)
    message("\nAlternative methods:")
    message("1. Use CELLxGENE web interface: https://cellxgene.cziscience.com/")
    message("2. Use cellxgenedp directly:")
    message("   library(cellxgenedp)")
    message("   datasets <- datasets()")
    message("   files <- files()")
    message("   files_download(file_id, destdir = 'output')")
    
    return(list(
      success = FALSE,
      error = e$message
    ))
  })
}
