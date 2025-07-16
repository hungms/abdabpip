#' Find hamming distance
#'
#' @param query The query data frame
#' @param reference The reference data frame
#' @param cols_to_match The columns to match
#' @param ncores The number of cores to use
#' @return A data frame with the hamming distance
#' @import future furrr stringdist alakazam dplyr stringr rlang magrittr data.table
#' @export
find_hamming_dist <- function(query, reference, cols_to_match, ncores){

    message(paste0("\nRunning hamming distance matching..."))

    # Setup future backend for parallelization
    future::plan(future::multisession, workers = ncores)

    # get CDR3 and gene columns to match
    cdr3_to_match <- cols_to_match[str_detect(names(cols_to_match), "CDR3")]
    genes_to_match <- cols_to_match[!str_detect(names(cols_to_match), "CDR3")]

    # Pre-filter query by gene matching if needed (do this once)
    if(length(genes_to_match) > 0){
        # Create a key for gene combinations to avoid repeated filtering
        query[, gene_key := do.call(paste, c(.SD, sep = "|")), .SDcols = names(genes_to_match)]
        reference[, gene_key := do.call(paste, c(.SD, sep = "|")), .SDcols = paste0("ref_", names(genes_to_match))]
        
        # Get unique gene combinations that exist in both query and reference
        valid_gene_keys <- intersect(query$gene_key, reference$gene_key)
        
        if(length(valid_gene_keys) == 0){
            message("No matching gene combinations found between query and reference")
            return(data.table())
        }
        
        # Filter both datasets to only include valid gene combinations
        query <- query[gene_key %in% valid_gene_keys]
        reference <- reference[gene_key %in% valid_gene_keys]
    }

    # Pre-calculate CDR3 lengths for query (do this once)
    for(col in names(cdr3_to_match)){
        query[, paste0(col, "_length") := nchar(get(col))]
    }

    # Using future + furrr with optimized approach
    output <- furrr::future_map_dfr(1:nrow(reference), function(i) {
        
        # Get reference row data once
        ref_row <- reference[i, ]
        
        # Filter query by gene matching (if needed) - much faster now
        if(length(genes_to_match) > 0){
            tmp <- query[gene_key == ref_row$gene_key]
        } else {
            tmp <- query
        }
        
        # Early exit if no matches
        if(nrow(tmp) == 0) return(data.table())
        
        # Add reference columns efficiently
        ref_cols <- paste0("ref_", names(cols_to_match))
        for(col in ref_cols){
            tmp[[col]] <- ref_row[[col]]
        }
        
        # Process CDR3 columns with length filtering
        for(col in names(cdr3_to_match)){
            ref_length_col <- paste0("ref_", col, "_length")
            query_length_col <- paste0(col, "_length")
            
            # Filter by matching CDR3 length
            tmp <- tmp[get(query_length_col) == ref_row[[ref_length_col]]]
            
            # Early exit if no matches
            if(nrow(tmp) == 0) return(data.table())
            
            # Get reference CDR3 sequence
            ref_cdr3 <- ref_row[[paste0("ref_", col)]]
            
            # Validate reference CDR3
            if(!str_detect(ref_cdr3, "^(?=(?:[^A-Z]*[A-Z]){5})[A-Z]+$")){
                tmp[, paste0(col, "_dist") := NA]
                tmp[, paste0(col, "_dist_method") := "hamming"]
                next
            }
            
            # Vectorized distance calculation for all matching sequences
            query_cdr3s <- tmp[[col]]
            valid_mask <- !is.na(query_cdr3s) & 
                         str_detect(query_cdr3s, "^(?=(?:[^A-Z]*[A-Z]){5})[A-Z]+$")
            
            if(any(valid_mask)){
                # Use vectorized stringdist for better performance
                distances <- stringdist::stringdist(ref_cdr3, query_cdr3s[valid_mask], method = "hamming") / nchar(ref_cdr3)
                tmp[valid_mask, paste0(col, "_dist") := distances]
                tmp[!valid_mask, paste0(col, "_dist") := NA]
            } else {
                tmp[, paste0(col, "_dist") := NA]
            }
            
            tmp[, paste0(col, "_dist_method") := "hamming"]
        }

        return(tmp)
    }, .options = furrr::furrr_options(seed = TRUE))

    # Clean up temporary columns
    if(length(genes_to_match) > 0){
        output[, gene_key := NULL]
    }

    return(output)
}

#' Find levenshtein distance
#'
#' @param query The query data frame
#' @param reference The reference data frame
#' @param cols_to_match The columns to match
#' @param ncores The number of cores to use
#' @return A data frame with the levenshtein distance
#' @import future furrr stringdist alakazam dplyr stringr rlang magrittr
#' @export
find_levenshtein_dist <- function(query, reference, cols_to_match, ncores){

    message(paste0("\nRunning levenshtein distance matching..."))

    # Setup future backend for parallelization
    future::plan(future::multisession, workers = ncores)

    # get CDR3 and gene columns to match
    cdr3_to_match <- cols_to_match[str_detect(names(cols_to_match), "CDR3")]
    genes_to_match <- cols_to_match[!str_detect(names(cols_to_match), "CDR3")]

    # Pre-filter query by gene matching if needed (do this once)
    if(length(genes_to_match) > 0){
        # Create a key for gene combinations to avoid repeated filtering
        query[, gene_key := do.call(paste, c(.SD, sep = "|")), .SDcols = names(genes_to_match)]
        reference[, gene_key := do.call(paste, c(.SD, sep = "|")), .SDcols = paste0("ref_", names(genes_to_match))]
        
        # Get unique gene combinations that exist in both query and reference
        valid_gene_keys <- intersect(query$gene_key, reference$gene_key)
        
        if(length(valid_gene_keys) == 0){
            message("No matching gene combinations found between query and reference")
            return(data.table())
        }
        
        # Filter both datasets to only include valid gene combinations
        query <- query[gene_key %in% valid_gene_keys]
        reference <- reference[gene_key %in% valid_gene_keys]
    }

    # Pre-calculate CDR3 lengths for query (do this once)
    for(col in names(cdr3_to_match)){
        query[, paste0(col, "_length") := nchar(get(col))]
    }

    # Using future + furrr
    output <- furrr::future_map_dfr(1:nrow(reference), function(i) {
        
        # add CDR3 sequence and length columns from reference data frame in query data frame
        tmp <- copy(query)
        for(col in colnames(reference)){
            tmp[[col]] <- reference[i, ][[col]]
        }
            
        # filter query for matching CDR3 sequence in heavy/light chains
        if(length(genes_to_match) > 0){
            for(col in names(genes_to_match)){
                tmp <- tmp[get(col) == get(paste0("ref_", col))]
            }
        }

        # levenshtein distance does not require matching CDR3 length
        for(col in names(cdr3_to_match)){

            # check if there are any matching sequences
            if(nrow(tmp) == 0) {next}

            # calculate levenshtein distance
            y <- reference[i, ][[paste0("ref_", col)]]
            tmp[, paste0(col, "_dist_method") := "levenshtein"]

            # query must have valid CDR3 sequence
            cond1 <- !str_detect(y, "^(?=(?:[^A-Z]*[A-Z]){5})[A-Z]+$")
            cond2 <- any(is.na(tmp[[col]]))
                
            if(any(cond1, cond2)){
                tmp[, paste0(col, "_dist") := NA]
            } else {
                tmp[, paste0(col, "_dist") := as.numeric(lapply(get(col), function(x){
                    stringdist::stringdist(y, x, method = 'lv') / max(nchar(y), nchar(x))
                }))]
            }
        }

        return(tmp)
    }, .options = furrr::furrr_options(seed = TRUE))

    return(output)
}

#' Find minimum distances
#'
#' @param output The output data frame
#' @param cols_to_match The columns to match
#' @param output_dir The output directory
#' @return A data frame with the minimum 
#' @import dplyr stringr magrittr
#' @export
find_min_distances <- function(output, cols_to_match, output_dir = NULL){

    message(paste0("Finding minimum CDR3 distance for each barcode..."))
    
    # Convert to data.table if not already
    if (!is.data.table(output)) {
        output <- as.data.table(output)
    }
    
    # remove rows with NA in ref_heavyCDR3
    output <- output[!is.na(ref_heavyCDR3)]

    # filter for rows where all dist_method values are the same
    dist_cols <- colnames(output)[str_detect(colnames(output), "_dist$")]
    dist_method_cols <- colnames(output)[str_detect(colnames(output), "_dist_method$")]
    
    # Check if all distance methods are the same for each row
    output[, all_dist_methods_same := uniqueN(.SD), .SDcols = dist_method_cols, by = 1:nrow(output)]
    output <- output[all_dist_methods_same == 1]
    output[, all_dist_methods_same := NULL]
    
    # check if there are any matching sequences
    if(nrow(output) > 0){

        # Calculate mean distance for each row
        output[, mean_dist := rowMeans(.SD, na.rm = TRUE), .SDcols = dist_cols]
        
        # Find minimum distance for each barcode
        output <- output[output[, .I[which.min(mean_dist)], by = barcodes]$V1]

        # found matching BCR sequences
        message(paste0("Found ", nrow(output), " public BCR sequences in QUERY..."))

        #write output to file
        if(!is.null(output_dir)){

            heavy_name <- names(cols_to_match)[str_detect(names(cols_to_match), "heavy")]
            heavy_name <- rev(sort(heavy_name))
            heavy_name <- paste0("heavy", paste0(gsub("heavy", "", heavy_name), collapse = ""))

            light_name <- names(cols_to_match)[str_detect(names(cols_to_match), "light")]
            light_name <- rev(sort(light_name))
            light_name <- paste0("light", paste0(gsub("light", "", light_name), collapse = ""))

            filename <- paste0(output_dir, "/publicBCR_by_", heavy_name, "_", light_name, ".csv")
            write.csv(output, filename, row.names = F)}

        # return output
        return(output)}

    # if no matching sequences, return NULL
    else{
        warning("\nNo public BCR sequences found")
        return(NULL)}
}