#' Detect Public BCR Sequences
#' 
#' @param query The query data frame
#' @param reference The reference data frame
#' @param heavyCDR3 The heavy chain CDR3 sequence
#' @param heavyV The heavy chain V gene
#' @param heavyJ The heavy chain J gene
#' @param lightCDR3 The light chain CDR3 sequence
#' @param lightV The light chain V gene
#' @param lightJ The light chain J gene
#' @param dist_method The distance method to use
#' @param ncores The number of cores to use
#' @param output_dir The output directory
#' @return A data frame with the matched CDR3 sequences
#' @import dplyr stringr magrittr data.table
#' @export

find_publicBCR <- function(
    query,
    reference,
    heavyCDR3,
    heavyV = NA,
    heavyJ = NA,
    lightCDR3 = NA,
    lightV = NA,
    lightJ = NA,
    dist_method = c("levenshtein", "hamming"),
    ncores = 1,
    output_dir = NULL){

    # Determine columns to match
    #========================================================
    cols_to_match <- lookup_cols(query, reference, heavyCDR3, heavyV, heavyJ, lightCDR3, lightV, lightJ)
    
    # Run preflight checks for query and reference
    #========================================================
    preflight_checks(ncores, dist_method, output_dir)
    query <- preflight_query(query, cols_to_match)
    reference <- preflight_reference(reference, cols_to_match)

    message(paste0("\nMatching ", nrow(query), " BCR sequences in QUERY against ", nrow(reference), " BCR sequences in REFERENCE..."))

    # Run Gene and CDR3 Matching
    #========================================================
    # create list to store outputs  
    output.list <- list()

    # Find hamming distance
    if("hamming" %in% dist_method){
        output.list[[length(output.list) + 1]] <- find_hamming_dist(query, reference, cols_to_match, ncores)}

    # Find levenshtein distance
    if("levenshtein" %in% dist_method){
        output.list[[length(output.list) + 1]] <- find_levenshtein_dist(query, reference, cols_to_match, ncores)}

    # combine outputs
    output <- bind_rows(output.list)

    # Find minimum distances
    #========================================================
    find_min_distances(output, cols_to_match, output_dir)
}
