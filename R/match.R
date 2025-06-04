#' Match CDR3 sequences between query and antigen collection
#' 
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

    # check ncores is valid
    stopifnot(ncores > 0 & is.numeric(ncores))

    # check antigen is valid
    stopifnot(antigen %in% c("Sars-CoV-2", "Tetanus", "Vaccinia", "Measles", "Mumps"))

    # check org is valid
    stopifnot(org %in% c("human", "mouse"))

    # check dist_method is valid
    stopifnot(all(dist_method %in% c("levenshtein", "hamming")))

    # check output_dir is valid
    stopifnot(dir.exists(output_dir))

    # check there are columns to match
    genes_to_match <- c(heavyV, heavyJ, lightV, lightJ)
    names(genes_to_match) <- c("heavyV", "heavyJ", "lightV", "lightJ")
    CDR3_to_match <- c(heavyCDR3, lightCDR3)
    names(CDR3_to_match) <- c("heavyCDR3", "lightCDR3")
    cols_to_match <- c(genes_to_match, CDR3_to_match)
    stopifnot(any(!is.na(cols_to_match)))

    # check columns to match are in query data frame
    cols_to_match <- na.omit(cols_to_match)
    stopifnot(all(cols_to_match %in% colnames(query)))

    # Logging
    #========================================================
    log_message()

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
        mutate(
            ref_VDJ_CDR3_length = nchar(VDJ_CDR3),
            ref_VJ_CDR3_length = nchar(VJ_CDR3))

    # Run VJ Match
    #========================================================
    # register parallel backend
    doParallel::registerDoParallel(ncores)

    # create list to store outputs
    output.list <- list()

    # if hamming is in dist_method, run matching
    if("hamming" %in% dist_method){

        # for each reference sequence
        output.list[[1]] <- foreach::foreach(i = 1:nrow(reference), .combine=rbind) %dopar% {

            # merge reference data frame to each query sequence by GENE columnes
            if(any(!is.na(genes_to_match))){
                tmp <- query %>%
                    semi_join(reference[i, ], by = c(names(genes_to_match)))}
                    
            # otherwise, add CDR3 sequence and length columns from reference data frame in query data frame
            else{
                tmp <- query
                for(col in c(names(CDR3_to_match), paste0(names(CDR3_to_match), "_length"))){
                    tmp[[col]] <- reference[i, ][[col]]}
            }

            # hamming distance requires matching CDR3 length
            if(any(str_detect(cols_to_match, "CDR3"))){
                for(col in cols_to_match[str_detect(cols_to_match, "CDR3")]){

                    # filter query for matching CDR3 length in heavy/light chains
                    tmp <- tmp %>% 
                        filter(!!sym(paste0(col, "_length")) == !!sym(paste0("ref_", col, "_length")))

                    # check if there are any matching sequences
                    if(nrow(tmp) == 0) {
                        warning("No matching BCR sequences found")
                        return(data.frame())}

                    # calculate hamming distance
                    y <- tmp[[paste0("ref_", col)]]
                    tmp[[paste0(col, "_hamming_dist")]] <- as.numeric(lapply(tmp[[col]], function(x){alakazam::seqDist(y, x, alakazam::getAAMatrix()) / nchar(y)}))}}
                
            return(tmp)}


    # if levenshtein is in dist_method, run matching
    if("levenshtein" %in% dist_method){

        # for each reference sequence
        output.list[[2]] <- foreach::foreach(i = 1:nrow(reference), .combine=rbind) %dopar% {

            # merge reference data frame to each query sequence by GENE columnes
            if(any(!is.na(genes_to_match))){
                tmp <- query %>%
                    semi_join(reference[i, ], by = c(names(genes_to_match)))}
                    
            # otherwise, add CDR3 sequence and length columns from reference data frame in query data frame
            else{
                tmp <- query
                for(col in c(names(CDR3_to_match), paste0(names(CDR3_to_match), "_length"))){
                    tmp[[col]] <- reference[i, ][[col]]}}

            #  levenshtein distance does not require matching CDR3 length
            if(any(str_detect(cols_to_match, "CDR3"))){
                for(col in cols_to_match[str_detect(cols_to_match, "CDR3")]){

                    # check if there are any matching sequences
                    if(nrow(tmp) == 0) {
                        warning("No matching BCR sequences found")
                        return(data.frame())}

                    # calculate hamming distance
                    y <- tmp[[paste0("ref_", col)]]
                    tmp[[paste0(col, "_levenshtein_dist")]] <- as.numeric(lapply(tmp$CDR3, function(x){stringdist::stringdist(cdr3, x, method = 'lv') / max(nchar(cdr3), nchar(x))}))}
                    }
            
            return(tmp)}
        }}

    # combine outputs
    output <- merge(output.list[[1]], output.list[[2]], by = "barcodes", all = T)

    # Format Output
    #========================================================
    # check if there are any matching sequences
    if(nrow(output) > 0){

        # get all distance columns
        dist_cols <- grepl("_dist$", colnames(output))

        # order output data frame 
        output <- output %>% 
            group_by(barcodes) %>% # group by barcode
            summarise(across(all_of(dist_cols), min, .names = "min_{.col}")) %>% # get minimum distance for each distance column
            return(output)

        # write output to file
        if(!is.null(output_dir)){
            filename <- paste0(output_dir, "/", antigen, "_", org, "_", "match_", paste0(cols_to_match, collapse = "_"), ".csv")
            write.csv(output, filename, row.names = F)}

        # return output
        return(output)}

    # if no matching sequences, return NULL
    else{
        warning("No matching BCR sequences found")
        return(NULL)}
}