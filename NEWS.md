# detectBCR 0.0.2

- Major refactor: split matching logic into modular functions in R/find_distances.R, R/find_publicBCR.R, R/preflight_checks.R, and R/lookup_cols.R
- All matching and preflight functions now use data.table for efficient processing
- Improved gene and CDR3 filtering logic for both query and reference
- Pre-filtering of query and reference by gene combinations for speed
- Vectorized hamming distance calculation for significant speedup
- Improved error handling and NA management in preflight and matching steps
- Removed redundant and batch/benchmark code for a cleaner codebase
- Updated vignette and README to reflect new function structure and usage

### v0.0.1
replace antigen and org arguments with 

### v0.0.0
