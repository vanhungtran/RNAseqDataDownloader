.onLoad <- function(libname, pkgname) {
  # Nothing needed here - all packages loaded via NAMESPACE
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("RNAseqDataDownloader loaded successfully. Use available_databases() to see supported sources.")
}
