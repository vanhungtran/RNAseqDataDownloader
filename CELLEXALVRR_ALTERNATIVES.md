# Replacement for Obsolete cellexalvrR Package

## Overview

The `cellexalvrR` package is **obsolete** and no longer maintained. This document provides modern alternatives for accessing Human Cell Atlas and single-cell data.

## Modern Alternative: cellxgenedp

**Package**: `cellxgenedp` (Bioconductor)  
**Purpose**: Access data from CZ CELLxGENE Discover portal  
**Status**: Active, maintained by Bioconductor  
**Data**: Human Cell Atlas, Mouse Cell Atlas, and other single-cell datasets

### Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("cellxgenedp")
```

### cellexalvrR vs cellxgenedp Function Mapping

| cellexalvrR (Obsolete) | cellxgenedp (Current) | Description |
|------------------------|----------------------|-------------|
| `getDatasets()` | `datasets()` | List available datasets |
| `getFiles()` | `files()` | List files for datasets |
| `downloadData()` | `files_download()` | Download dataset files |
| `getCollections()` | `collections()` | Browse data collections |
| `getCells()` | Via dataset metadata | Get cell count info |

## Updated Function: download_hca()

The `download_hca()` function has been completely rewritten to use `cellxgenedp`:

### New Features

1. **Browse datasets** with filters
2. **Download H5AD** (AnnData) or RDS (Seurat) formats
3. **Filter by organism**, tissue, or disease
4. **Full metadata** access

### Usage Examples

#### 1. Browse Available Datasets

```r
library(RNAseqDataDownloader)

# Show available datasets (filtered)
result <- download_hca(
  organism = "Homo sapiens",
  tissue = "lung"
)

# View datasets
head(result$available_datasets)
```

#### 2. Download a Specific Dataset

```r
# Download by dataset ID
result <- download_hca(
  dataset_id = "b9fc3d70-5a72-4479-a046-c2cc1ab19efc",
  output_dir = "hca_data",
  file_format = "h5ad"
)

# File location
print(result$file_path)
```

#### 3. Load and Analyze Downloaded Data

```r
# For H5AD files - use anndata or Seurat
library(Seurat)
library(SeuratDisk)

# Convert H5AD to Seurat
Convert(result$file_path, dest = "h5seurat", overwrite = TRUE)
seurat_obj <- LoadH5Seurat(gsub("\\.h5ad$", ".h5seurat", result$file_path))

# Or download directly as RDS
result <- download_hca(
  dataset_id = "your-dataset-id",
  file_format = "rds"
)
seurat_obj <- readRDS(result$file_path)
```

## Direct cellxgenedp Usage

For more control, use `cellxgenedp` directly:

### Browse Datasets

```r
library(cellxgenedp)

# Get all datasets
all_datasets <- datasets()
head(all_datasets)

# Filter datasets
human_lung <- all_datasets |>
  dplyr::filter(
    organism == "Homo sapiens",
    grepl("lung", tissue, ignore.case = TRUE)
  )

# View dataset details
View(human_lung[, c("dataset_id", "dataset_title", "cell_count", "disease")])
```

### Download Files

```r
# Get files for a dataset
all_files <- files()
dataset_files <- all_files |>
  dplyr::filter(dataset_id == "your-dataset-id")

# Download a specific file
local_path <- files_download(
  file_id = dataset_files$file_id[1],
  destdir = "output_dir",
  dry.run = FALSE  # Set TRUE to see what would be downloaded
)
```

### Get Collections

```r
# Browse curated collections
collections_data <- collections()
head(collections_data[, c("collection_id", "name", "description")])

# Get datasets in a collection
collection_datasets <- collections() |>
  dplyr::filter(collection_id == "your-collection-id")
```

## Other Single-Cell Data Sources

### 1. Single Cell Expression Atlas (EBI)

Still available via `download_sc_atlas()`:

```r
# Download from EBI Single Cell Expression Atlas
result <- download_sc_atlas(
  accession = "E-MTAB-5061",
  output_dir = "sc_data",
  file_type = "loom"
)
```

### 2. GEO Single-Cell Data

Many single-cell studies are in GEO:

```r
# Download single-cell data from GEO
result <- download_geo(
  "GSE123456",
  data_type = "raw",
  filter_files = "\\.h5$|\\.mtx$"  # Filter for single-cell formats
)
```

### 3. HCA Data Portal (Manual)

For official Human Cell Atlas data:

- **Website**: https://data.humancellatlas.org/
- **CLI Tool**: https://github.com/HumanCellAtlas/dcp-cli
- **API**: https://data.humancellatlas.org/apis

```bash
# Install HCA CLI
pip install hca

# List projects
hca dss search --query "organ.text:lung"

# Download files
hca dss download --manifest manifest.json
```

## Migration Guide

### Old Code (cellexalvrR)

```r
# OLD - Don't use!
library(cellexalvrR)

datasets <- getDatasets()
files <- getFiles(dataset_id)
downloadData(file_id, output_dir)
```

### New Code (cellxgenedp)

```r
# NEW - Use this!
library(cellxgenedp)

datasets <- datasets()
files_list <- files()
dataset_files <- files_list[files_list$dataset_id == dataset_id, ]
local_path <- files_download(dataset_files$file_id[1], destdir = output_dir)
```

### Using RNAseqDataDownloader Wrapper

```r
# Even simpler - use our wrapper
library(RNAseqDataDownloader)

# Browse datasets
available <- download_hca(organism = "Homo sapiens", tissue = "brain")

# Download specific dataset
result <- download_hca(
  dataset_id = "your-dataset-id",
  output_dir = "hca_data"
)
```

## Advantages of cellxgenedp

1. ✅ **Actively maintained** by Bioconductor
2. ✅ **Large collection** - 1000+ datasets
3. ✅ **Well documented** with vignettes
4. ✅ **Multiple formats** - H5AD, RDS, CSV
5. ✅ **Rich metadata** - cell types, markers, protocols
6. ✅ **Fast downloads** with caching
7. ✅ **Quality controlled** data
8. ✅ **Regular updates** with new datasets

## Additional Resources

### cellxgenedp Documentation

- **Bioconductor**: https://bioconductor.org/packages/cellxgenedp
- **Vignette**: https://bioconductor.org/packages/devel/bioc/vignettes/cellxgenedp/inst/doc/using_cellxgenedp.html
- **GitHub**: https://github.com/mtmorgan/cellxgenedp

### CZ CELLxGENE Portal

- **Web Interface**: https://cellxgene.cziscience.com/
- **Documentation**: https://chanzuckerberg.github.io/cellxgene/
- **API Docs**: https://api.cellxgene.cziscience.com/

### Human Cell Atlas

- **Portal**: https://data.humancellatlas.org/
- **Documentation**: https://www.humancellatlas.org/learn/
- **Data Coordination Platform**: https://data.humancellatlas.org/apis

## Example Workflow

Complete workflow using the new approach:

```r
library(RNAseqDataDownloader)
library(Seurat)
library(dplyr)

# 1. Browse datasets
available <- download_hca(
  organism = "Homo sapiens",
  tissue = "pancreas",
  disease = "type 2 diabetes"
)

# 2. Select a dataset
dataset_id <- available$available_datasets$dataset_id[1]
print(available$available_datasets[1, c("dataset_title", "cell_count")])

# 3. Download
result <- download_hca(
  dataset_id = dataset_id,
  output_dir = "my_analysis/data",
  file_format = "rds"  # Get as Seurat-compatible RDS
)

# 4. Load and analyze
seurat_obj <- readRDS(result$file_path)

# 5. QC and analysis
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# 6. Visualization
DimPlot(seurat_obj, group.by = "cell_type")
```

## Summary

| Feature | cellexalvrR | cellxgenedp |
|---------|-------------|-------------|
| Status | ❌ Obsolete | ✅ Active |
| Maintenance | ❌ None | ✅ Bioconductor |
| Datasets | Limited | 1000+ |
| Formats | Limited | H5AD, RDS, CSV |
| Documentation | ❌ Outdated | ✅ Current |
| Community | ❌ Inactive | ✅ Active |

**Recommendation**: Use `cellxgenedp` for all single-cell data downloads from CELLxGENE and Human Cell Atlas.
