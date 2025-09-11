.onLoad <- function(libname, pkgname) {
  # Check if required Bioconductor packages are installed
  required_bioc_pkgs <- c("GEOquery", "TCGAbiolinks", "recount3", "SRAdb", "ArrayExpress")

  for(pkg in required_bioc_pkgs) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing required Bioconductor package: ", pkg)
      BiocManager::install(pkg)
    }
  }

  # Load the packages
  suppressPackageStartupMessages({
    library(GEOquery, quietly = TRUE)
    library(TCGAbiolinks, quietly = TRUE)
    library(recount3, quietly = TRUE)
    library(SRAdb, quietly = TRUE)
    library(ArrayExpress, quietly = TRUE)
  })

  message("RNAseqDataDownloader loaded successfully. Use available_databases() to see supported sources.")
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to RNAseqDataDownloader!")
}
