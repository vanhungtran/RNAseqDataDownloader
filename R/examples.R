#' Example usage of RNAseqDataDownloader
#' @export
example_usage <- function() {
  # Set a larger timeout for big downloads
  set_timeout(600)

  # Set default download directory
  set_download_dir("my_rnaseq_data")

  # Download GEO data
  geo_matrix <- download_geo("GSE12345", data_type = "matrix")

  # Download SRA data
  sra_files <- download_sra("SRR1234567", method = "fasterq")

  # Download TCGA data
  tcga_data <- download_tcga("TCGA-BRCA")

  # Download GTEx data
  gtex_counts <- download_gtex(data_type = "counts")

  message("Example downloads completed!")
}
