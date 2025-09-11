#' Download data from SRA
#' @param accession SRA accession number (e.g., "SRR1234567")
#' @param output_dir Directory to save downloaded data
#' @param method Download method: "fasterq" (default) or "prefetch"
#' @return Path to downloaded files
#' @export
download_sra <- function(accession, output_dir = "sra_data", method = "fasterq") {
  check_and_create_dir(output_dir)

  # Check if SRA Toolkit is installed
  if (Sys.which("prefetch") == "" || Sys.which("fasterq-dump") == "") {
    stop("SRA Toolkit is not installed or not in PATH. Please install it from https://github.com/ncbi/sra-tools")
  }

  tryCatch({
    if (method == "prefetch") {
      message("Using prefetch to download ", accession)
      system(paste("prefetch", accession, "-O", output_dir))
      sra_file <- file.path(output_dir, accession, paste0(accession, ".sra"))

      # Convert to FASTQ
      message("Converting SRA to FASTQ")
      system(paste("fasterq-dump", sra_file, "-O", file.path(output_dir, accession)))
      return(file.path(output_dir, accession))
    } else if (method == "fasterq") {
      message("Using fasterq-dump to download ", accession)
      system(paste("fasterq-dump", accession, "-O", output_dir))
      return(file.path(output_dir, paste0(accession, ".fastq")))
    } else {
      stop("method must be 'fasterq' or 'prefetch'")
    }
  }, error = function(e) {
    stop("Failed to download SRA data: ", e$message)
  })
}

#' Get SRA metadata
#' @param accession SRA accession number
#' @return Data frame with metadata
#' @export
get_sra_metadata <- function(accession) {
  # Use ENA API as SRAdb is being deprecated
  ena_url <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
                    accession, "&result=read_run&fields=all&format=json")

  response <- httr::GET(ena_url)
  if (httr::status_code(response) != 200) {
    stop("Failed to retrieve metadata for ", accession)
  }

  metadata <- jsonlite::fromJSON(httr::content(response, "text"))
  return(metadata)
}
