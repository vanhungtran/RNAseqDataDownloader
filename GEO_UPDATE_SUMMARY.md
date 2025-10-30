# GEO Download Function Update Summary

## Overview
Updated the `download_geo()` function with latest best practices from the GEOquery GitHub repository (seandavi/GEOquery, versions 2.99.0+, 2024-2025).

## Major Updates

### 1. **RNA-seq Quantification Support**
- Added new `data_type = "rnaseq"` option to download NCBI-processed RNA-seq count data
- Uses GEOquery's `getRNASeqData()` function for human and mouse datasets
- Returns SummarizedExperiment objects with raw counts and metadata
- Automatically falls back to matrix download if RNA-seq data unavailable

### 2. **Annotation GPL Support**
- Added `AnnotGPL` parameter (default: FALSE)
- When TRUE, downloads GPL files with updated Entrez Gene mappings
- Provides more current gene annotations than submitter GPLs
- Available for GPLs referenced by GDS datasets

### 3. **Enhanced GSEMatrix Performance**
- Uses GSEMatrix=TRUE by default (10-100x faster than SOFT format)
- Improved handling of multi-platform GSE entries
- Better platform filtering and selection
- More efficient memory usage

### 4. **Improved Download Management**
- Added `timeout` parameter (default: 300 seconds)
- Configurable timeout for large file downloads
- Better retry logic with exponential backoff
- Enhanced error handling and fallback mechanisms

### 5. **Better File Handling**
- Updated `filter_regex` parameter for getGEOSuppFiles
- Support for tar archives extraction
- More robust file format detection
- Saves both text and RDS formats for easy re-loading

### 6. **Enhanced Metadata Extraction**
- More comprehensive platform information
- Better handling of characteristics parsing
- Added `parseCharacteristics` parameter
- Improved sample filtering support

## New Features

### RNA-seq Data Download
```r
# Download RNA-seq quantification data
result <- download_geo("GSE164073", data_type = "rnaseq")

# Returns SummarizedExperiment with:
# - Raw counts matrix
# - Gene annotations
# - Sample metadata
# - Genome build information
```

### Annotation GPL Usage
```r
# Use updated gene annotations
result <- download_geo("GSE12345", AnnotGPL = TRUE)
```

### Advanced File Filtering
```r
# Download only specific file types
result <- download_geo("GSE12345", 
                      data_type = "raw",
                      filter_files = "\\.txt\\.gz$")
```

## Updated Parameters

### New Parameters
- `AnnotGPL` - Use annotation GPL with updated mappings
- `timeout` - Download timeout in seconds
- `parseCharacteristics` - Parse characteristics fields

### Modified Parameters
- `data_type` - Added "rnaseq" option, removed "miniml"
- `filter_files` - Now uses GEOquery's native filtering

### Enhanced Return Values
- `platform_info` - Detailed platform information
- `download_time` - Timestamp of download
- `success` - Boolean success indicator

## Implementation Details

### Based on GEOquery GitHub
- Repository: seandavi/GEOquery
- Latest features from version 2.99.0+ (2024)
- Uses httr2 for HTTP requests
- Improved error handling patterns
- Better memory management

### Key Changes
1. **GSEMatrix First**: Prioritizes fast GSEMatrix format
2. **RNA-seq Detection**: Automatically detects available RNA-seq data
3. **Robust Retries**: Exponential backoff for network issues
4. **Platform Filtering**: Better multi-platform GSE handling
5. **Archive Extraction**: Enhanced support for tar, tar.gz, gz, zip

## Backwards Compatibility

### Maintained
- All original parameters still work
- Existing code should run without changes
- Default behavior unchanged (matrix download)

### Removed
- `miniml` data type (rarely used, not in latest GEOquery)

## Usage Examples

### Basic Download (Fastest)
```r
result <- download_geo("GSE12345")
```

### RNA-seq Data
```r
result <- download_geo("GSE164073", data_type = "rnaseq")
```

### Complete Download
```r
result <- download_geo("GSE12345", 
                      data_type = "both",
                      AnnotGPL = TRUE,
                      timeout = 600)
```

### Platform-Specific
```r
result <- download_geo("GSE12345",
                      platform = "GPL570",
                      data_type = "matrix")
```

### Sample Filtering
```r
result <- download_geo("GSE12345",
                      samples = c("GSM123", "GSM124"),
                      data_type = "both")
```

## Performance Improvements

1. **GSEMatrix Parsing**: 10-100x faster than SOFT format
2. **Parallel Downloads**: Better handling of large files
3. **Memory Efficiency**: Improved data structure handling
4. **Caching**: Respects local cache to avoid re-downloads

## Error Handling

- Better error messages with specific guidance
- Automatic fallback mechanisms
- Exponential backoff for network issues
- Detailed verbose output for troubleshooting

## Testing Recommendations

Test with various accession types:
```r
# GSE with RNA-seq
download_geo("GSE164073", data_type = "rnaseq")

# Multi-platform GSE
download_geo("GSE12345", data_type = "matrix")

# Platform annotation
download_geo("GPL570", AnnotGPL = TRUE)

# Dataset
download_geo("GDS507")

# With supplementary files
download_geo("GSE12345", data_type = "both", 
            filter_files = "\\.cel\\.gz$")
```

## References

- GEOquery GitHub: https://github.com/seandavi/GEOquery
- GEO RNA-seq Documentation: https://ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html
- GEOquery Vignette: https://bioconductor.org/packages/GEOquery
- NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/

## Future Enhancements

Potential additions based on GEOquery development:
1. Single-cell data support (h5ad, mtx formats)
2. Direct search integration
3. Batch download optimization
4. SummarizedExperiment conversion utilities
5. Integration with other Bioconductor packages

---

**Updated**: October 30, 2025
**GEOquery Version Reference**: 2.99.0+
**Maintainer**: RNAseqDataDownloader package
