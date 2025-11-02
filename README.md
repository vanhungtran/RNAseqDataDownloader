
- [Tutorial for R package
  ‘RNAseqDataDownloader’](#tutorial-for-r-package-rnaseqdatadownloader)
  - [Introduction](#introduction)
- [🧬 RNAseqDataDownloader](#dna-rnaseqdatadownloader)
  - [*Download ARNSeq Public Data in
    R*](#download-arnseq-public-data-in-r)
- [📖 Overview](#open_book-overview)
  - [✅ Key Features](#white_check_mark-key-features)
- [🚀 Installation](#rocket-installation)
- [🧪 Quick Examples](#test_tube-quick-examples)
- [📚 Core Functions](#books-core-functions)
- [📊 Examples](#bar_chart-examples)
- [🤝 Contributing](#handshake-contributing)
- [📜 License](#scroll-license)
- [❓ Need Help?](#question-need-help)
- [Reference](#reference)
- [🙏 Acknowledgements](#pray-acknowledgements)

<!-- README.md is auto-generated from README.Rmd -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- README.md is auto-generated from README.Rmd -->

# Tutorial for R package ‘RNAseqDataDownloader’

Lucas TRAN 10/09/2025

## Introduction

*RNAseqDataDownloader*

<div align="center">

# 🧬 RNAseqDataDownloader

### *Download ARNSeq Public Data in R*

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/yourusername/RNAseqDataDownloader/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourusername/RNAseqDataDownloader/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/RNAseqDataDownloader)](https://cran.r-project.org/package=RNAseqDataDownloader)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<br>

> *“Classify, validate, and explore drugs with WHO ATC codes — all in
> R.”*

</div>

------------------------------------------------------------------------

## 📖 Overview

**`RNAseqDataDownloader`** is an R package that simplifies working with
**ATC codes** — the World Health Organization’s global classification
system for medicinal products. Whether you’re analyzing prescription
databases, mapping drug classes, validating inputs, or linking drug
names to codes, `RNAseqDataDownloader` provides intuitive, fast, and
reliable tools.

### ✅ Key Features

- ✔ High-level wrappers to download GEO, TCGA, GTEx, SRA, HCA, and single-cell atlas datasets
- ✔ Unified metadata extraction with phenotype parsing, platform details, and SRA links
- ✔ Batch download tooling with retry logic, logging, and matrix post-processing helpers
- ✔ Configurable download directory, timeouts, and archive extraction for reproducible pipelines
- ✔ Export utilities to convert GEO objects into tidy data frames, matrices, or SummarizedExperiment
- ✔ Vignettes and examples covering end-to-end RNA-seq acquisition workflows

Ideal for **transcriptomics**, **single-cell analysis**, **biomedical data integration**, and
**reproducible RNA-seq pipelines**.

------------------------------------------------------------------------

## 🚀 Installation

The code of *RNAseqDataDownloader* is freely available at
<https://github.com/vanhungtran/RNAseqDataDownloader>.

The following commands can be used to install this R package, and an R
version \>= 4.2.3 is required.

    library(devtools)
    install_github("vanhungtran/RNAseqDataDownloader")

Load the package:

``` r
library(RNAseqDataDownloader)
```

## 🧪 Quick Examples

## 📚 Core Functions

| Function                    | Description                                                                 |
|-----------------------------|-----------------------------------------------------------------------------|
| `available_databases()`     | List all public repositories currently supported by the downloader          |
| `download_geo()`            | Download GEO series/platform/sample datasets (matrix, raw, RNA-seq, SOFT)   |
| `batch_download_geo()`      | Batch download multiple GEO accessions with retry handling and summaries    |
| `extract_geo_metadata()`    | Parse rich metadata (phenotype, platforms, SRA links) from GEO objects      |
| `download_tcga()`           | Retrieve TCGA RNA-seq data and metadata via the GDC API                     |
| `download_gtex()`           | Download GTEx expression matrices, metadata, and sample annotations         |
| `download_sc_atlas()`       | Access curated single-cell atlases with unified metadata                    |
| `download_hca()`            | Fetch Human Cell Atlas data with optional filtering and caching             |
| `download_sra()`            | Download SRA runs linked to studies or GEO accessions                       |
| `set_download_dir()`        | Configure a persistent default directory for all dataset downloads          |
| `load_geo_batch_matrices()` | Load every GEO series matrix file in a folder into a named list of matrices |
| `combine_geo_matrices()`    | Merge matrices into a long-format table and add a `GEOcode` identifier      |
| `load_geo_batch_pheno()`    | Read and clean GEO phenotype tables across batch downloads with unified IDs |

------------------------------------------------------------------------

## 📊 Examples

Let’s explore how

``` r
```

------------------------------------------------------------------------

## 🤝 Contributing

We welcome contributions! Please see our
[CONTRIBUTING.md](CONTRIBUTING.md) guide for how to:

- Report bugs
- Suggest features
- Submit pull requests
- Improve documentation

------------------------------------------------------------------------

## 📜 License

MIT © 2025 \[Lucas VHH Tran\]

> Permission is hereby granted, free of charge, to any person obtaining
> a copy of this software…

See [LICENSE.md](LICENSE.md) for full text.

------------------------------------------------------------------------

------------------------------------------------------------------------

## ❓ Need Help?

Open an issue on
[GitHub](https://github.com/vanhungtran/RNAseqDataDownloader/issues) or
email \[<tranhungydhcm@gmail.com>\].

------------------------------------------------------------------------

<div align="center">

<br> <em>Developed with RNAseq and 🧬 for the R and health data science
community.</em> <br><br>

</div>

------------------------------------------------------------------------

<div align="center">

<img src="man/figures/logo.png" width="220">

</div>

## Reference

## 🙏 Acknowledgements

- GEOquery, Biobase, and GEOmetadb teams for maintaining the GEO access
  ecosystem
- National Cancer Institute GDC, GTEx Consortium, and Human Cell Atlas
  for providing open access RNA-seq resources
- Contributors to Bioconductor packages (`SummarizedExperiment`,
  `GenomicDataCommons`, `SingleCellExperiment`, `BiocParallel`) used in
  this package
- Open-source community maintaining supporting tools (`httr2`, `jsonlite`,
  `arrow`, `data.table`, `cli`)

<!---
&#10;
4. Commit both `README.Rmd` and `README.md` to GitHub.
&#10;> 💡 Tip: Add `README.md` to `.Rbuildignore` if you don’t want it in the built package (optional).
&#10;
## 🖌️ Customize Further
- Replace `yourusername` with your GitHub username
- Replace `[Your Name or Organization]` and email
- Add your hex sticker image under `man/figures/logo.png` and uncomment the image line if desired
- Update example code to match your actual function names and outputs
&#10;
&#10;
&#10;🎁 Bonus: Generate a Hex Sticker in R
&#10;If you want to create a hex sticker, install `hexSticker` and run:
&#10;
library(ggplot2)
library(hexSticker)
&#10;# Just a centered "A" for ATC
p <- ggplot() + 
  annotate("text", x = 1, y = 1, label = "A", size = 20, fontface = "bold") +
  xlim(0.5, 1.5) + ylim(0.5, 1.5) +
  theme_void()
&#10;sticker(
  subplot = p,
  package = "RNAseqDataDownloader",
  p_size = 10,
  s_x = 1,
  s_y = 0.8,
  s_width = 1.3,
  filename = "man/figures/logo.png",
  h_fill = "#2a9d8f",
  h_color = "#264653",
&#10;)
 &#10;-->
