# detectBCR Dependencies

This document outlines all dependencies for the detectBCR package, their purposes, and installation instructions.

## Required Dependencies (Imports)

These packages are automatically installed when you install detectBCR:

### Core Data Manipulation
- **dplyr**: Data manipulation and transformation functions
- **magrittr**: Pipe operator (`%>%`) for chaining operations
- **tibble**: Modern data frame implementation with better printing
- **stringr**: String manipulation and pattern matching

### Parallel Processing
- **future**: Unified interface for parallel processing
- **furrr**: Parallel processing with purrr-like syntax

### Sequence Analysis
- **alakazam**: Immunoglobulin sequence analysis and distance calculations
- **stringdist**: String distance calculations (Levenshtein, Hamming, etc.)

## Optional Dependencies (Suggests)

These packages provide performance enhancements and are automatically used if available:

### Performance Optimization
- **data.table**: High-performance data manipulation (5-20x speedup for post-processing)
- **gpuR**: GPU acceleration for distance calculations (5-20x speedup for large datasets)

### Documentation and Testing
- **knitr**: Dynamic report generation
- **rmarkdown**: R Markdown document generation
- **testthat**: Unit testing framework

## Installation

### Basic Installation
```r
install.packages("devtools")
devtools::install_github("hungms/detectBCR", dependencies = TRUE)
```

### With Performance Enhancements
```r
# Install optional performance packages
install.packages(c("data.table", "gpuR"))

# Install detectBCR
devtools::install_github("hungms/detectBCR", dependencies = TRUE)
```

### GPU Requirements (for gpuR)
- **CUDA Toolkit**: Version 8.0 or higher
- **GPU**: CUDA-compatible NVIDIA GPU
- **R**: Version 3.5.0 or higher

## Dependency Usage in detectBCR

### Core Functions
- **Data filtering**: `dplyr`, `stringr`
- **Distance calculations**: `alakazam`, `stringdist`
- **Parallel processing**: `future`, `furrr`

### Performance Features
- **GPU acceleration**: `gpuR` (automatic detection)
- **Fast data processing**: `data.table` (automatic fallback to dplyr)

### Documentation
- **Vignettes**: `knitr`, `rmarkdown`
- **Testing**: `testthat`

## Version Compatibility

The package is tested with:
- R >= 3.5.0
- dplyr >= 1.0.0
- future >= 1.20.0
- furrr >= 0.2.0
- alakazam >= 1.0.0
- stringdist >= 0.9.0

## Troubleshooting

### GPU Issues
If GPU acceleration fails:
1. Check CUDA installation: `nvidia-smi`
2. Verify gpuR installation: `library(gpuR)`
3. Package will automatically fallback to CPU

### Performance Issues
For slow performance:
1. Install `data.table`: `install.packages("data.table")`
2. Install `gpuR` (if GPU available): `install.packages("gpuR")`
3. Use more cores: `ncores = parallel::detectCores() * 0.75`

### Memory Issues
For large datasets:
1. Use fewer cores to reduce memory usage
2. Process data in chunks
3. Use `data.table` for memory-efficient operations 