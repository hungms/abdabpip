pkgs <- c(
    "dplyr", "magrittr", "stringr", "future", "furrr", "alakazam", "stringdist", "tibble")

for(x in pkgs){
    usethis::use_package(x, type = "depends")} #, type = "depends"