#' Match CDR3 sequences between query and antigen collection
#' 
#' @param query The query data frame
#' @param antigen The antigen to match
#' @param org The organism to match
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
#' @export

match_CDR3 <- function(
    query,
    antigen,
    org,
    heavyCDR3,
    heavyV = NA,
    heavyJ = NA,
    lightCDR3 = NA,
    lightV = NA,
    lightJ = NA,
    dist_method = c("levenshtein", "hamming"),
    ncores = 1,
    output_dir = NULL){
    

    # Validate Inputs
    #========================================================
    # check query is a data frame
    stopifnot(is.data.frame(query))

    # check if there are any invalid CDR3 sequences in heavyCDR3
    ## remove rows with less than 5 amino acids in heavyCDR3
    rows_to_remove <- !str_detect(query[[heavyCDR3]], "^(?=(?:[^A-Z]*[A-Z]){5})[A-Z]+$")
    if(any(rows_to_remove)){
        message(paste0("\nRemoving ", sum(rows_to_remove), " rows with invalid CDR3 sequences in ", heavyCDR3, " in QUERY data frame..."))
        query <- query[!rows_to_remove, ]}

    # check ncores is valid
    stopifnot(ncores > 0 & is.numeric(ncores))

    # check antigen is valid
    stopifnot(antigen %in% c("Sars-CoV-2", "Tetanus", "Vaccinia", "Measles", "Mumps"))

    # check org is valid
    stopifnot(org %in% c("human", "mouse"))

    # check dist_method is valid
    stopifnot(all(dist_method %in% c("levenshtein", "hamming")))

    # check output_dir is valid
    if(!is.null(output_dir)){
        stopifnot(dir.exists(output_dir))}

    # store columns to match
    genes_to_match <- c(heavyV, heavyJ, lightV, lightJ)
    names(genes_to_match) <- c("heavyV", "heavyJ", "lightV", "lightJ")
    CDR3_to_match <- c(heavyCDR3, lightCDR3)
    names(CDR3_to_match) <- c("heavyCDR3", "lightCDR3")

    # check columns to match are in query data frame
    genes_to_match <- genes_to_match[!is.na(genes_to_match)]
    CDR3_to_match <- CDR3_to_match[!is.na(CDR3_to_match)]
    cols_to_match <- c(genes_to_match, CDR3_to_match)

    stopifnot(all(cols_to_match %in% colnames(query)))

    # Logging
    #========================================================
    #log_message()

    # Format Data
    #========================================================
    # change colnames in query data frame
    for(col in cols_to_match){
        colnames(query)[which(colnames(query) == col)] <- names(cols_to_match)[cols_to_match == col]}

    # select columns to match from query and reference data frames
    query <- query %>%
        rownames_to_column("barcodes") %>%
        dplyr::select(barcodes, names(cols_to_match))

    # get collection
    collection <- get_reference(antigen = antigen)

    # add CDR3 length columns in query data frame
    if(any(!is.na(CDR3_to_match))){
        for(i in names(CDR3_to_match)){
            query <- query %>%
                mutate(!!paste0(i, "_length") := nchar(!!sym(i)))}}

    # add CDR3 length columns in reference data frame
    reference <- collection %>%
        filter(ref_org == org) %>%
        filter(ref_heavyCDR3 != "") %>%
        filter(!is.na(ref_heavyCDR3)) %>%
        mutate(
            ref_heavyCDR3_length = nchar(ref_heavyCDR3),
            ref_lightCDR3_length = nchar(ref_lightCDR3))

    message(paste0("\nMatching ", nrow(query), " BCR sequences in QUERY against ", nrow(reference), " BCR sequences in ", org, " REFERENCE..."))
    message(paste0("\nMatching columns: ", paste0(names(cols_to_match), collapse = ", "), "..."))
    message(paste0("\nSearching for ", antigen, "-specific BCR sequences..."))

    # Run VJ Match
    #========================================================
    # Setup future backend for parallelization
    future::plan(future::multisession, workers = ncores)
    
    # create list to store outputs
    output.list <- list()

    # if hamming is in dist_method, run matching
    if("hamming" %in% dist_method){
        message(paste0("\nRunning hamming distance matching..."))

        # Using future + furrr
        output.list[[length(output.list) + 1]] <- furrr::future_map_dfr(1:nrow(reference), function(i) {
            # add CDR3 sequence and length columns from reference data frame in query data frame
            tmp <- query
            for(col in colnames(reference)){
                tmp[[col]] <- reference[i, ][[col]]
            }

            # filter query for matching CDR3 sequence in heavy/light chains
            if(length(genes_to_match) > 0){
                for(col in names(genes_to_match)){
                    tmp <- tmp %>%
                        filter(!!sym(paste0(col)) == !!sym(paste0("ref_", col)))
                }
            }
            
            # hamming distance requires matching CDR3 length
            for(col in names(cols_to_match)[str_detect(names(cols_to_match), "CDR3")]){

                # filter query for matching CDR3 length in heavy/light chains
                tmp <- tmp %>% 
                    filter(!!sym(paste0(col, "_length")) == !!sym(paste0("ref_", col, "_length")))

                # check if there are any matching sequences
                if(nrow(tmp) == 0) {
                    next
                }

                # calculate hamming distance
                y <- reference[i, ][[paste0("ref_", col)]]
                tmp[[paste0(col, "_dist_method")]] <- "hamming"

                # query must have valid CDR3 sequence
                cond1 <- !str_detect(y, "^(?=(?:[^A-Z]*[A-Z]){5})[A-Z]+$")
                cond2 <- any(is.na(tmp[[col]]))

                if(any(cond1, cond2)){
                    tmp[[paste0(col, "_dist")]] <- NA
                } else {
                    tmp[[paste0(col, "_dist")]] <- as.numeric(lapply(tmp[[col]], function(x){
                        alakazam::seqDist(y, x, alakazam::getAAMatrix()) / nchar(y)
                    }))
                }
            }

            return(tmp)
        }, .options = furrr::furrr_options(seed = TRUE))
    }

    # if levenshtein is in dist_method, run matching
    if("levenshtein" %in% dist_method){
        message(paste0("\nRunning levenshtein distance matching..."))
        
        # Using future + furrr
        output.list[[length(output.list) + 1]] <- furrr::future_map_dfr(1:nrow(reference), function(i) {
            # add CDR3 sequence and length columns from reference data frame in query data frame
            tmp <- query
            for(col in colnames(reference)){
                tmp[[col]] <- reference[i, ][[col]]
            }
            
            # filter query for matching CDR3 sequence in heavy/light chains
            if(length(genes_to_match) > 0){
                for(col in names(genes_to_match)){
                    tmp <- tmp %>%
                        filter(!!sym(paste0(col)) == !!sym(paste0("ref_", col)))
                }
            }

            #  levenshtein distance does not require matching CDR3 length
            for(col in names(cols_to_match)[str_detect(names(cols_to_match), "CDR3")]){

                # check if there are any matching sequences
                if(nrow(tmp) == 0) {
                    next
                }

                # calculate levenshtein distance
                y <- reference[i, ][[paste0("ref_", col)]]
                tmp[[paste0(col, "_dist_method")]] <- "levenshtein"

                # query must have valid CDR3 sequence
                cond1 <- !str_detect(y, "^(?=(?:[^A-Z]*[A-Z]){5})[A-Z]+$")
                cond2 <- any(is.na(tmp[[col]]))
                
                if(any(cond1, cond2)){
                    tmp[[paste0(col, "_dist")]] <- NA
                } else {
                    tmp[[paste0(col, "_dist")]] <- as.numeric(lapply(tmp[[col]], function(x){
                        stringdist::stringdist(y, x, method = 'lv') / max(nchar(y), nchar(x))
                    }))
                }
            }
            
            return(tmp)
        }, .options = furrr::furrr_options(seed = TRUE))
    }

    # combine outputs
    output <- bind_rows(output.list)

    # get all distance columns
    dist_cols <- colnames(output)[str_detect(colnames(output), "_dist$")]
    dist_method_cols <- colnames(output)[str_detect(colnames(output), "_dist_method$")]

    # combine outputs
    output <- output %>%
        # remove rows with NA in ref_heavyCDR3
        filter(!is.na(ref_heavyCDR3)) 
    
    # Filter for rows where all dist_method values are the same
    output <- output %>%
        rowwise() %>%
        mutate(all_dist_methods_same = n_distinct(c_across(all_of(dist_method_cols)))) %>%
        filter(all_dist_methods_same == 1) %>%
        select(-all_dist_methods_same) %>%
        ungroup()

    # Format Output
    #========================================================
    # check if there are any matching sequences
    if(nrow(output) > 0){
        message(paste0("Finding minimum CDR3 distance for each barcode..."))

        # order output data frame 
        output <- output %>% 
            rowwise() %>%
            mutate(sum_dist = mean(c_across(all_of(dist_cols)), na.rm = T)) %>%
            ungroup() %>%
            group_by(barcodes) %>% # group by barcode
            slice_min(n = 1, order_by = sum_dist) %>% # get the barcode with the minimum sum of distances
            ungroup()

        # found matching BCR sequences
        message(paste0("Found ", nrow(output), " BCR sequences in QUERY that has a match in the REFERENCE..."))

        #write output to file
        if(!is.null(output_dir)){
            filename <- paste0(output_dir, "/", antigen, "_", org, "_", "match_", paste0(names(cols_to_match), collapse = "_"), ".csv")
            write.csv(output, filename, row.names = F)}

        # return output
        return(output)}

    # if no matching sequences, return NULL
    else{
        warning("\nNo matching BCR sequences found")
        return(NULL)}
}