#' Download GTEx data
#' @param output_dir Directory to save downloaded data
#' @param data_type Type of data to download: "counts", "tpm", or "clinical"
#' @return Path to downloaded files
#' @export
download_gtex <- function(output_dir = "gtex_data", data_type = "counts") {
  check_and_create_dir(output_dir)

  tryCatch({
    # Use recount3 to access GTEx data
    proj_info <- recount3::available_projects()
    
    # Fix: Add explicit column checks
    project_type <- NULL
    project_home <- NULL
    gtex_proj <- subset(proj_info, proj_info$project_type == "data_sources" & proj_info$project_home == "gtex")

    if (nrow(gtex_proj) == 0) {
      stop("GTEx data not available in recount3")
    }

    message("Downloading GTEx ", data_type, " data")

    # Create a URL for the specific data type
    base_url <- "https://storage.googleapis.com/gtex_analysis_v8/"

    urls <- list(
      counts = paste0(base_url, "rna_seq_data/gene_reads/gene_reads_2017-06-05_v8.gct.gz"),
      tpm = paste0(base_url, "rna_seq_data/gene_tpm/gene_tpm_2017-06-05_v8.gct.gz"),
      clinical = paste0(base_url, "reference/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt")
    )

    if (!data_type %in% names(urls)) {
      stop("data_type must be one of: ", paste(names(urls), collapse = ", "))
    }

    dest_file <- file.path(output_dir, basename(urls[[data_type]]))
    download.file(urls[[data_type]], destfile = dest_file)

    message("GTEx data downloaded successfully")
    return(dest_file)
  }, error = function(e) {
    stop("Failed to download GTEx data: ", e$message)
  })
}
