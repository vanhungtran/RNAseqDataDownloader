library(testthat)

test_that("download_tcga dry-run works (mocked)", {
  skip_if_not_installed("mockery")
  skip_if_not_installed("SummarizedExperiment")

  # Minimal mocks to avoid network calls
  dummy_query <- list(files = data.frame(cases = I(list()), file_id = "dummy", stringsAsFactors = FALSE))
  qmock <- function(...) dummy_query
  dmock <- function(...) invisible(TRUE)
  pmock <- function(...) {
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = matrix(1, 1)),
      colData = S4Vectors::DataFrame(sample = 1)
    )
    se
  }

  # Stub external TCGAbiolinks functions used inside download_tcga
  mockery::stub(RNAseqDataDownloader::download_tcga, 'TCGAbiolinks::GDCquery', qmock)
  mockery::stub(RNAseqDataDownloader::download_tcga, 'TCGAbiolinks::GDCdownload', dmock)
  mockery::stub(RNAseqDataDownloader::download_tcga, 'TCGAbiolinks::GDCprepare', pmock)

  res <- NULL
  expect_error_free(res <- RNAseqDataDownloader::download_tcga(
    project = "TCGA-FAKE",
    prepare_data = FALSE,
    download_clinical = FALSE
  ))

  expect_true(is.list(res) || is.null(res))
})
