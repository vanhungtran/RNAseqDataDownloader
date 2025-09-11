# usethis::use_description(
#   list(
#     Package = "RNAseqDataDownloader",
#     Title = "Unified Interface for Downloading RNA-seq Data from Public Databases",
#     Version = "0.1.0",
#     Authors@R = person("Lucas", "VHH Tran", email = "tranhungzdhcmgmail.com", role = c("aut", "cre")),
#     Description = "Provides a consistent interface for downloading RNA-seq data from various public databases including GEO, SRA, TCGA, GTEx, and single-cell atlases.",
#     License = "MIT",
#     Depends = "R (>= 4.0.0)",
#     Imports = "GEOquery, TCGAbiolinks, recount3, SRAdb, ArrayExpress, cellexalvrR, Seurat, httr, jsonlite, data.table, BiocManager",
#     Suggests = "testthat, knitr, rmarkdown",
#     VignetteBuilder = "knitr"
#   )
# )


# Install usethis if not already installed
if (!requireNamespace("usethis", quietly = TRUE)) {
  install.packages("usethis")
}

# Create a new package
usethis::create_package("~/RNAseqDataDownloader")

# Set the working directory to the package
setwd("~/RNAseqDataDownloader")

# Initialize git repository
usethis::use_git()
usethis::use_github() # If you want to use GitHub

use_mit_license("Lucas VHH Tran")


# Set up package dependencies
usethis::use_package("GEOquery")
usethis::use_package("TCGAbiolinks")
usethis::use_package("recount3")
usethis::use_package("SRAdb")
usethis::use_package("ArrayExpress")
usethis::use_package("httr")
usethis::use_package("jsonlite")
usethis::use_package("data.table")
usethis::use_package("BiocManager")
usethis::use_package("Seurat")
usethis::use_package("cellexalvrR")





usethis::use_r("geo")
usethis::use_r("sra")
usethis::use_r("tcga")
usethis::use_r("gtex")
usethis::use_r("sc_atlas")
usethis::use_r("utils")
usethis::use_r("zzz")


# Add Bioconductor repositories
usethis::use_description_field("BioCViews", "AnnotationData Software")

# We also need to set the repository for Bioconductor
options(repos = BiocManager::repositories())





















# Developing an R Package for RNA-seq Data Download using usethis

creating an R package for downloading RNA-seq data from various databases using the `usethis` package.
This approach follows modern R package development practices.

## Step 1: Setting up the package structure

First, let's create the basic package structure:

  ```r
# Install usethis if not already installed
if (!requireNamespace("usethis", quietly = TRUE)) {
  install.packages("usethis")
}

# Create a new package
usethis::create_package("~/RNAseqDataDownloader")

# Set the working directory to the package
setwd("~/RNAseqDataDownloader")

# Initialize git repository
usethis::use_git()
usethis::use_github() # If you want to use GitHub

# Add license
usethis::use_mit_license("Your Name")

# Set up package dependencies
usethis::use_package("GEOquery")
usethis::use_package("TCGAbiolinks")
usethis::use_package("recount3")
usethis::use_package("SRAdb")
usethis::use_package("ArrayExpress")
usethis::use_package("httr")
usethis::use_package("jsonlite")
usethis::use_package("data.table")
usethis::use_package("BiocManager")
usethis::use_package("Seurat")
usethis::use_package("cellexalvrR")

# For packages that are not on CRAN (like some Bioconductor packages)
# Add to DESCRIPTION manually:
# Imports:
#     GEOquery,
#     TCGAbiolinks,
#     recount3,
#     SRAdb,
#     ArrayExpress,
#     cellexalvrR,
#     Seurat,
#     httr,
#     jsonlite,
#     data.table,
#     BiocManager
```

## Step 2: Creating the R functions

Now let's create the R functions using `usethis`:

```r
# Create utility functions
usethis::use_r("utils")
```

In the newly created `R/utils.R` file:

```r
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
```

Create the GEO functions:

  ```r
usethis::use_r("geo")
```

In `R/geo.R`:

  ```r
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
```

Create the SRA functions:

  ```r
usethis::use_r("sra")
```

In `R/sra.R`:

  ```r
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
```

Create the TCGA functions:

  ```r
usethis::use_r("tcga")
```

In `R/tcga.R`:

  ```r
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
```

Create the GTEx functions:

  ```r
usethis::use_r("gtex")
```

In `R/gtex.R`:

  ```r
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
    gtex_proj <- subset(proj_info, project_type == "data_sources" & project_home == "gtex")

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
```

Create the single-cell functions:

  ```r
usethis::use_r("sc_atlas")
```

In `R/sc_atlas.R`:

  ```r
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
```

Create the configuration functions:

  ```r
usethis::use_r("config")
```

In `R/config.R`:

  ```r
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
```

## Step 3: Package initialization

Create the package initialization file:

  ```r
usethis::use_r("zzz")
```

In `R/zzz.R`:

  ```r
.onLoad <- function(libname, pkgname) {
  # Check if required Bioconductor packages are installed
  required_bioc_pkgs <- c("GEOquery", "TCGAbiolinks", "recount3", "SRAdb", "ArrayExpress")

  for(pkg in required_bioc_pkgs) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing required Bioconductor package: ", pkg)
      BiocManager::install(pkg)
    }
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to RNAseqDataDownloader! Use available_databases() to see supported sources.")
}
```

## Step 4: Documentation and testing

Create a README file:

  ```r
usethis::use_readme_rmd()
```

Create vignettes:

  ```r
usethis::use_vignette("introduction")
usethis::use_vignette("geo_download")
usethis::use_vignette("tcga_download")
```

Set up testing:

  ```r
usethis::use_testthat()
usethis::use_test("utils")
usethis::use_test("geo")
# Add tests for other modules as needed
```

## Step 5: Build and check the package

```r
# Build the package
devtools::document() # Generate documentation
devtools::check() # Check for errors
devtools::build() # Build the package

# Install the package
devtools::install()
```

## Step 6: Example usage

Create an example script:

  ```r
usethis::use_r("examples")
```

In `R/examples.R`:

  ```r
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
```

## Final steps

1. Add more tests to ensure package reliability
2. Consider adding a function to check SRA Toolkit installation
3. Add error handling for various network conditions
4. Consider adding progress bars for large downloads
5. Add more examples to vignettes

This package structure provides a solid foundation for downloading RNA-seq data from various sources. The `usethis` package helps maintain a consistent and professional package structure.
