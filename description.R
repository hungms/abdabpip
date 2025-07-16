pkgs <- c(
    "dplyr", "magrittr", "stringr", "future", "furrr", "stringdist", "tibble", "data.table")

for(x in pkgs){
    usethis::use_package(x, type = "depends")}