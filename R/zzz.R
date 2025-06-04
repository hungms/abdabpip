.onLoad <- function(libname, pkgname) {
  ### load default packages
  packages <- c(
    "dplyr", "magrittr", "stringr", "foreach", "doParallel", "alakazam", "stringdist")
  invisible(lapply(packages, library, character.only = TRUE))

  ### start up settings
  options(dplyr.summarise.inform = FALSE)
  set.seed(123)
  errorMessage <- NULL
}