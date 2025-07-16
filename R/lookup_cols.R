#' Lookup columns to match
#'
#' @param query The query data frame
#' @param reference The reference data frame
#' @param heavyCDR3 The heavy chain CDR3 sequence
#' @param heavyV The heavy chain V gene
#' @param heavyJ The heavy chain J gene
#' @param lightCDR3 The light chain CDR3 sequence
#' @param lightV The light chain V gene
#' @param lightJ The light chain J gene
#' @return A list of columns to match
#' @import dplyr magrittr stringr data.table
#' @export
lookup_cols <- function(
    query,
    reference,
    heavyCDR3,
    heavyV = NA,
    heavyJ = NA,
    lightCDR3 = NA,
    lightV = NA,
    lightJ = NA){

    # store gene columns to match
    genes_to_match <- c(heavyV, heavyJ, lightV, lightJ)
    names(genes_to_match) <- c("heavyV", "heavyJ", "lightV", "lightJ")
    genes_to_match <- genes_to_match[!is.na(genes_to_match)]

    # store CDR3 columns to match
    CDR3_to_match <- c(heavyCDR3, lightCDR3)
    names(CDR3_to_match) <- c("heavyCDR3", "lightCDR3")
    CDR3_to_match <- CDR3_to_match[!is.na(CDR3_to_match)]

    # check if all columns to match are in query and reference
    cols_to_match <- c(genes_to_match, CDR3_to_match)
    stopifnot(all(cols_to_match %in% colnames(query)))
    stopifnot(all(paste0("ref_", names(cols_to_match)) %in% colnames(reference)))

    return(cols_to_match)
    }