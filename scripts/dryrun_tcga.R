#!/usr/bin/env Rscript
# Safe dry-run script for download_tcga()
# Usage (PowerShell):
#   Rscript scripts/dryrun_tcga.R                # runs default dry-run (project TCGA-FAKE)
#   Rscript scripts/dryrun_tcga.R TCGA-BRCA      # runs with project override (still a dry-run by default)

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args) >= 1) args[[1]] else "TCGA-FAKE"

# Non-destructive defaults for safety
prepare_data <- FALSE
download_clinical <- FALSE
download_biospecimen <- FALSE
max_retries <- 1
verbose <- TRUE

options(stringsAsFactors = FALSE)

# Ensure devtools is available to load package from source
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

# Load package from current working directory (the repo root)
devtools::load_all(".")

# Run a safe dry-run of download_tcga()
res <- tryCatch(
  {
    download_tcga(
      project = project,
      prepare_data = prepare_data,
      download_clinical = download_clinical,
      download_biospecimen = download_biospecimen,
      max_retries = max_retries,
      verbose = verbose
    )
  },
  error = function(e) list(error = conditionMessage(e))
)

# Save result for inspection
out_file <- file.path(".", paste0("dryrun_tcga_result_", gsub("[^A-Za-z0-9_-]", "_", project), ".rds"))
saveRDS(res, out_file)
cat("Dry-run complete. Result saved to:", normalizePath(out_file), "\n")

# Helpful note for real runs (commented):
# To perform an actual download, set prepare_data = TRUE and/or
# download_clinical = TRUE, download_biospecimen = TRUE. Example:
# download_tcga(project = "TCGA-BRCA", prepare_data = TRUE, download_clinical = TRUE)
