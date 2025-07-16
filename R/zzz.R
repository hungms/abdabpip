.onLoad <- function(libname, pkgname) {
  ### load default packages
  packages <- c(
    "dplyr", "magrittr", "stringr", "future", "furrr", "stringdist", "tibble", "data.table")
  invisible(lapply(packages, library, character.only = TRUE))

  ### start up settings
  options(dplyr.summarise.inform = FALSE)
  set.seed(123)
}