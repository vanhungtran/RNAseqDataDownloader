# Create Complete R Package "mestools"
# This script creates a fully-featured R package with all components

# Load required packages
required_packages <- c("usethis", "devtools", "desc", "pkgdown", "knitr", "rmarkdown")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}))

library(usethis)
library(devtools)

# Set package name
pkg_name <- "mestools"

# Create package structure
message("Creating package structure for '", pkg_name, "'...")

# Create basic package
create_package(path = pkg_name, rstudio = TRUE, open = FALSE)

# Set working directory to package
setwd(pkg_name)

# 1. BASIC PACKAGE SETUP
message("1. Setting up basic package structure...")

# Package description with all fields
use_description(
  list(
    Package = pkg_name,
    Title = "A Collection of Useful Tools for Data Analysis and Package Development",
    Description = "Provides utility functions for data manipulation, package development, deployment automation, and common analytical tasks. Includes GitHub deployment automation and quality of life tools for R developers.",
    `Authors@R` = c(
      person("Your", "Name", email = "your.email@example.com", role = c("aut", "cre"))
    ),
    Version = "0.1.0",
    License = "MIT + file LICENSE",
    URL = paste0("https://github.com/vanhungtran/", pkg_name),
    BugReports = paste0("https://github.com/vanhungtran/", pkg_name, "/issues"),
    Depends = "R (>= 3.5.0)",
    Imports = "utils, stats, tools, devtools, usethis, gert, here, desc, pkgdown, knitr, rmarkdown, testthat, covr",
    Suggests = "rmarkdown, knitr, spelling, pkgload, curl, gh",
    VignetteBuilder = "knitr",
    Encoding = "UTF-8"
  )
)

# 2. LICENSE
message("2. Adding MIT license...")
use_mit_license("Your Name")

# 3. DIRECTORY STRUCTURE
message("3. Creating directory structure...")

# Create standard directories
dir.create("R", showWarnings = FALSE, recursive = TRUE)
dir.create("man", showWarnings = FALSE, recursive = TRUE)
dir.create("vignettes", showWarnings = FALSE, recursive = TRUE)
dir.create("inst", showWarnings = FALSE, recursive = TRUE)
dir.create("tests/testthat", showWarnings = FALSE, recursive = TRUE)

# Create inst directories
dir.create("inst/doc", showWarnings = FALSE, recursive = TRUE)

# 4. CORE R FUNCTIONS
message("4. Creating core R functions...")

# Main package functions file
writeLines(
  '# Package: mestools
# Collection of utility functions

NULL

# Helper function for NULL coalescing
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

#' Validate GitHub Repository
#'
#' Checks if a GitHub repository exists and is accessible.
#'
#' @param repo_url GitHub repository URL
#' @return Logical indicating if repository exists
#' @export
#' @examples
#' \\dontrun{
#' validate_github_repo("https://github.com/r-lib/usethis.git")
#' }
validate_github_repo <- function(repo_url) {
  if (!requireNamespace("gh", quietly = TRUE)) {
    stop("Please install the \\'gh\\' package: install.packages(\\'gh\\')")
  }
  
  # Extract owner and repo name from URL
  repo_parts <- strsplit(repo_url, "/")[[1]]
  owner <- repo_parts[length(repo_parts)-1]
  repo_name <- gsub("\\\\.git$", "", repo_parts[length(repo_parts)])
  
  message("Checking GitHub repository: ", repo_url)
  
  tryCatch({
    repo_info <- gh::gh(paste0("/repos/", owner, "/", repo_name))
    message("✅ Repository exists: ", repo_info$html_url)
    message("   Description: ", repo_info$description %||% "No description")
    message("   Visibility: ", ifelse(repo_info$private, "Private", "Public"))
    return(TRUE)
  }, error = function(e) {
    if (grepl("404", e$message)) {
      message("❌ Repository not found: ", repo_url)
    } else {
      message("❌ Error checking repository: ", e$message)
    }
    return(FALSE)
  })
}

#' Deploy R Package to GitHub
#'
#' Complete workflow for rebuilding, testing, and deploying an R package to GitHub.
#'
#' @param repo_url GitHub repository URL (optional)
#' @param commit_message Commit message (optional)
#' @return Logical indicating success
#' @export
#' @examples
#' \\dontrun{
#' deploy_package()
#' }
deploy_package <- function(repo_url = NULL, commit_message = NULL) {
  # Load required packages
  required_packages <- c("devtools", "usethis", "gert", "here")
  invisible(lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }))
  
  package_name <- basename(here::here())
  
  if (is.null(repo_url)) {
    repo_url <- paste0("https://github.com/vanhungtran/", package_name, ".git")
  }
  
  if (is.null(commit_message)) {
    commit_message <- paste("Daily update for", package_name, "-", Sys.Date())
  }
  
  tryCatch({
    message("🚀 Starting deployment for: ", package_name)
    message("📦 Target repository: ", repo_url)
    
    # Validate GitHub repo first
    if (!validate_github_repo(repo_url)) {
      stop("GitHub repository validation failed. Deployment aborted.")
    }
    
    # Package rebuilding steps
    message("1. 📝 Cleaning and rebuilding documentation...")
    unlink("man", recursive = TRUE)
    dir.create("man")
    devtools::document()
    
    message("2. 📚 Rebuilding vignettes...")
    devtools::build_vignettes()
    
    message("3. 🧪 Running tests...")
    test_results <- devtools::test()
    if (any(as.data.frame(test_results)$failed > 0) {
      stop("Tests failed. Deployment aborted.")
    }
    
    message("4. ✅ Running package check...")
    check_results <- devtools::check(error_on = "never")
    if (length(check_results$errors) > 0 || length(check_results$warnings) > 0) {
      warning("Package check completed with warnings or errors - review manually")
    }
    
    # Git operations
    message("5. 🔗 Setting up git remote...")
    usethis::use_git_remote("origin", repo_url, overwrite = TRUE)
    
    message("6. 💾 Committing changes...")
    
    # Stage all changes
    gert::git_add(".")
    
    # Check if there are changes to commit
    status <- gert::git_status()
    if (nrow(status) == 0) {
      message("No changes to commit.")
      return(invisible(FALSE))
    }
    
    # Commit changes
    gert::git_commit(commit_message)
    
    message("7. 🚀 Pushing to GitHub...")
    gert::git_push(remote = "origin", set_upstream = TRUE)
    
    message("🎉 Deployment completed successfully!")
    return(invisible(TRUE))
    
  }, error = function(e) {
    message("❌ Deployment failed: ", e$message)
    return(invisible(FALSE))
  })
}

#' Quick Data Summary
#'
#' Provides a quick summary of a data frame.
#'
#' @param data Data frame to summarize
#' @return Summary information
#' @export
#' @examples
#' quick_summary(mtcars)
quick_summary <- function(data) {
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  list(
    dimensions = dim(data),
    column_names = names(data),
    column_types = sapply(data, class),
    missing_values = sapply(data, function(x) sum(is.na(x))),
    memory_size = format(utils::object.size(data), units = "auto")
  )
}

#' Safe File Reader
#'
#' Safely reads various file formats with error handling.
#'
#' @param file_path Path to the file
#' @param type File type (auto-detected if NULL)
#' @return File contents
#' @export
#' @examples
#' \\dontrun{
#' read_file_safe("data.csv")
#' }
read_file_safe <- function(file_path, type = NULL) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  if (is.null(type)) {
    type <- tools::file_ext(file_path)
  }
  
  tryCatch({
    switch(tolower(type),
           csv = read.csv(file_path),
           tsv = read.delim(file_path),
           rds = readRDS(file_path),
           txt = readLines(file_path),
           stop("Unsupported file type: ", type)
    )
  }, error = function(e) {
    stop("Error reading file: ", e$message)
  })
}
', "R/mestools.R")

# Utility functions file
writeLines(
'# Utility functions for mestools package

#' Check Package Dependencies
#'
#' Checks if all package dependencies are installed.
#'
#' @param packages Vector of package names
#' @return Logical indicating if all packages are available
#' @export
check_dependencies <- function(packages) {
  available <- sapply(packages, requireNamespace, quietly = TRUE)
  if (!all(available)) {
    missing <- packages[!available]
    message("Missing packages: ", paste(missing, collapse = ", "))
    return(FALSE)
  }
  TRUE
}

#' Create Project Directory Structure
#'
#' Creates a standardized directory structure for projects.
#'
#' @param path Base path for project
#' @param directories Vector of directory names to create
#' @export
create_project_structure <- function(path = ".", directories = c("data", "scripts", "output", "docs")) {
  for (dir in directories) {
    dir_path <- file.path(path, dir)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      message("Created: ", dir_path)
    }
  }
  message("Project structure created successfully!")
}

#' Batch Apply Function
#'
#' Applies a function to multiple objects.
#'
#' @param objects List of objects
#' @param fun Function to apply
#' @param ... Additional arguments to fun
#' @return List of results
#' @export
batch_apply <- function(objects, fun, ...) {
  lapply(objects, function(obj) {
    tryCatch({
      fun(obj, ...)
    }, error = function(e) {
      warning("Error processing object: ", e$message)
      NULL
    })
  })
}
', "R/utils.R")

# 5. TESTS
message("5. Setting up test suite...")

use_testthat()

# Test file for main functions
writeLines(
'test_that("NULL coalescing works", {
  expect_equal(NULL %||% "default", "default")
  expect_equal("value" %||% "default", "value")
})

test_that("quick_summary returns correct structure", {
  result <- quick_summary(mtcars)
  expect_named(result, c("dimensions", "column_names", "column_types", 
                         "missing_values", "memory_size"))
  expect_equal(result$dimensions, c(32, 11))
})

test_that("check_dependencies works", {
  # Should return TRUE for base packages
  expect_true(check_dependencies(c("utils", "stats")))
})
', "tests/testthat/test-mestools.R")

# 6. VIGNETTES
message("6. Creating vignettes...")

# Main vignette
writeLines(
'---
  title: "Getting Started with mestools"
author: "Your Name"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\\VignetteIndexEntry{Getting Started}
%\\VignetteEngine{knitr::rmarkdown}
%\\VignetteEncoding{UTF-8}
---
  
  # Introduction to mestools
  
  The `mestools` package provides a collection of utility functions for data analysis and package development.

## Key Features

### Package Deployment

```{r, eval=FALSE}
library(mestools)

# Deploy your package to GitHub
deploy_package()
