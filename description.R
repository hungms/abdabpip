pkgs <- c(
    "dplyr", "magrittr", "stringr", "foreach", "doParallel", "alakazam", "stringdist", "tibble")

for(x in pkgs){
    usethis::use_package(x, type = "depends")} #, type = "depends"