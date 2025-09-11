#' Download data from TCGA
#' @param project TCGA project ID (e.g., "TCGA-BRCA")
#' @param output_dir Directory to save downloaded data
#' @param data_category Data category to download (default: "Transcriptome Profiling")
#' @param data_type Data type to download (default: "Gene Expression Quantification")
#' @param workflow_type Workflow type (default: "HTSeq - Counts")
#' @return Path to downloaded files
#' @export
download_tcga <- function(project, output_dir = "tcga_data",
                          data_category = "Transcriptome Profiling",
                          data_type = "Gene Expression Quantification",
                          workflow_type = "HTSeq - Counts") {
  check_and_create_dir(output_dir)

  tryCatch({
    message("Querying TCGA data for project: ", project)

    # Build query
    query <- TCGAbiolinks::GDCquery(
      project = project,
      data.category = data_category,
      data.type = data_type,
      workflow.type = workflow_type
    )

    # Download data
    TCGAbiolinks::GDCdownload(query, directory = output_dir)

    # Prepare data
    data <- TCGAbiolinks::GDCprepare(query, directory = output_dir)

    # Save data as RDS for later use
    saveRDS(data, file.path(output_dir, paste0(project, "_data.rds")))

    message("TCGA data downloaded and prepared successfully")
    return(file.path(output_dir, paste0(project, "_data.rds")))
  }, error = function(e) {
    stop("Failed to download TCGA data: ", e$message)
  })
}
