Dry-run scripts for RNAseqDataDownloader

scripts/dryrun_tcga.R
- Purpose: run a safe, non-destructive dry-run of `download_tcga()`.
- Usage (PowerShell):

```powershell
cd 'C:\Users\tranh\OneDrive\RNAseqDataDownloader'
Rscript scripts/dryrun_tcga.R                # dry-run with default project (TCGA-FAKE)
Rscript scripts/dryrun_tcga.R TCGA-BRCA      # dry-run for a specific project
```

- Notes:
  - The script loads the package source via `devtools::load_all('.')`.
  - By default the script sets `prepare_data = FALSE` and disables clinical/biospecimen downloads to avoid network and large file operations.
  - To perform an actual download, run interactively and call `download_tcga(..., prepare_data = TRUE)` or edit the script accordingly.
