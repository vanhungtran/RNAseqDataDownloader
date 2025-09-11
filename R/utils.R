#' Check if output directory exists and create if not
#' @param output_dir Path to output directory
#' @return Logical indicating success
check_and_create_dir <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    message("Creating output directory: ", output_dir)
    return(dir.create(output_dir, recursive = TRUE))
  }
  return(TRUE)
}

#' List available databases supported by the package
#' @return Character vector of database names
#' @export
available_databases <- function() {
  c("GEO", "SRA", "TCGA", "GTEx", "SingleCellAtlas", "HCA", "ArrayExpress", "Recount3")
}

#' Set download timeout
#' @param seconds Timeout in seconds
#' @export
set_timeout <- function(seconds = 300) {
  options(timeout = seconds)
  message("Download timeout set to ", seconds, " seconds")
}



#' Check if output directory exists and create if not
#' @param output_dir Path to output directory
#' @return Logical indicating success
#' @keywords internal
check_and_create_dir <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    message("Creating output directory: ", output_dir)
    return(dir.create(output_dir, recursive = TRUE))
  }
  return(TRUE)
}

#' List available databases supported by the package
#' @return Character vector of database names
#' @export
available_databases <- function() {
  c("GEO", "SRA", "TCGA", "GTEx", "SingleCellAtlas", "HCA", "ArrayExpress", "Recount3")
}

#' Set download timeout
#' @param seconds Timeout in seconds
#' @export
set_timeout <- function(seconds = 300) {
  options(timeout = seconds)
  message("Download timeout set to ", seconds, " seconds")
}
