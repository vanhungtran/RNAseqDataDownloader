#' Download data from GEO
#' @param accession GEO accession number (e.g., "GSE12345")
#' @param output_dir Directory to save downloaded data
#' @param data_type Type of data to download: "matrix" (expression matrix) or "raw" (raw data)
#' @param ... Additional parameters passed to getGEO or getGEOSuppFiles
#' @return Path to downloaded files
#' @export
download_geo <- function(accession, output_dir = "geo_data", data_type = "matrix", ...) {
  check_and_create_dir(output_dir)

  tryCatch({
    if (data_type == "matrix") {
      message("Downloading GEO series matrix for ", accession)
      gse <- GEOquery::getGEO(accession, destdir = output_dir, ...)
      return(file.path(output_dir, paste0(accession, "_series_matrix.txt.gz")))
    } else if (data_type == "raw") {
      message("Downloading GEO supplementary files for ", accession)
      GEOquery::getGEOSuppFiles(accession, baseDir = output_dir, ...)
      return(file.path(output_dir, accession))
    } else {
      stop("data_type must be 'matrix' or 'raw'")
    }
  }, error = function(e) {
    stop("Failed to download GEO data: ", e$message)
  })
}

#' Extract metadata from GEO object
#' @param geo_object GEO object from getGEO
#' @return Data frame with metadata
#' @export
extract_geo_metadata <- function(geo_object) {
  if (!inherits(geo_object, "list") && !inherits(geo_object, "ExpressionSet")) {
    stop("Input must be a GEO object from getGEO()")
  }

  # Handle both list and single ExpressionSet
  if (inherits(geo_object, "list")) {
    pheno_data <- Biobase::pData(geo_object[[1]])
  } else {
    pheno_data <- Biobase::pData(geo_object)
  }

  return(pheno_data)
}
