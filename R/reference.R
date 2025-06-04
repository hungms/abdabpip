#' Get reference data frame
#' 
#' @param antigen The antigen to get the reference data frame for
#' @return A data frame with the reference data
#' @export

get_reference <- function(antigen){
    # check antigen is valid
    stopifnot(antigen %in% c("Sars-CoV-2", "Tetanus", "Vaccinia", "Measles", "Mumps"))
    reference <- read.csv(system.file("extdata", paste0(antigen, ".csv"), package = "detectBCR"), header = T, sep = ",")
    return(reference)
}
