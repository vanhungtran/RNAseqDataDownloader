#' Set default download directory
#' @param path Path to default download directory
#' @export
set_download_dir <- function(path = "rnaseq_data") {
  options(RNAseqDataDownloader.download_dir = path)
  check_and_create_dir(path)
  message("Default download directory set to: ", path)
}

#' Get default download directory
#' @return Path to default download directory
#' @export
get_download_dir <- function() {
  dir <- getOption("RNAseqDataDownloader.download_dir", "rnaseq_data")
  if (!dir.exists(dir)) {
    check_and_create_dir(dir)
  }
  return(dir)
}
