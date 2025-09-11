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

#' Download data from Human Cell Atlas
#' @param project_id HCA project ID
#' @param output_dir Directory to save downloaded data
#' @param file_type Type of file to download: "h5ad" or "loom"
#' @return Path to downloaded files
#' @export
download_hca <- function(project_id, output_dir = "hca_data", file_type = "h5ad") {
  check_and_create_dir(output_dir)

  tryCatch({
    # Use the HCA API or DCP CLI would be better, but for simplicity we'll use a direct approach
    # Note: This is a simplified example - actual HCA download would be more complex

    message("HCA download requires authentication and is complex.")
    message("Please use the HCA Data Browser or HCA CLI for downloading HCA data.")
    message("See: https://data.humancellatlas.org/")

    return(NULL)
  }, error = function(e) {
    stop("Failed to download HCA data: ", e$message)
  })
}
