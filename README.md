# AbDabList
Database & Pipeline to find matching CDR3 sequences from scVDJseq output (from Dandelion)


## Dependencies
```ruby
cran <- c("tidyverse", "foreach", "doParallel", "alakazam", "stringdist")
install.packages(cran)
```

## Find Matching CDR3

```ruby
library(tidyverse)
library(foreach)
library(doParallel)
library(alakazam)
library(stringdist)

match_aa <- function(query, reference, method = "hamming", match = c("v_call_VDJ_main", "j_call_VDJ_main"), ncores = 1){

  stopifnot(method %in% c("hamming", "levenshtein"))
  stopifnot(all(match %in% c("v_call_VDJ_main", "j_call_VDJ_main", "junction_aa_VDJ")))
  stopifnot(all(c(match, "junction_aa_VDJ") %in% colnames(query)))
  stopifnot(all(c(match, "junction_aa_VDJ") %in% colnames(reference)))

  query <- query %>%
    dplyr::select(c(match, "junction_aa_VDJ"))
  reference <- reference %>%
    distinct(!!!syms(c(match, "junction_aa_VDJ")))

  query$junction_aa_VDJ_length <- nchar(query$junction_aa_VDJ)
  reference$junction_aa_VDJ_length <- nchar(reference$junction_aa_VDJ)

  doParallel::registerDoParallel(ncores)

  output <- foreach::foreach(i = 1:nrow(reference), .combine=rbind) %dopar% {

    if(method == "hamming") {
      tmp <- query %>% semi_join(reference[i, ], by = c(match, "junction_aa_VDJ_length")) %>% rownames_to_column("cell_barcode")
      if(nrow(tmp) == 0) {
        return(data.frame())}
      cdr3 <- reference[i, ]$junction_aa_VDJ
      tmp$dist <- as.numeric(lapply(tmp$junction_aa_VDJ, function(x){alakazam::seqDist(cdr3, x, alakazam::getAAMatrix()) / nchar(cdr3)}))
      tmp$refseq <- cdr3
      return(tmp)}

    if(method == "levenshtein"){
      if(length(match) > 0) {
        tmp <- query %>% semi_join(reference[i, ], by = match) %>% rownames_to_column("cell_barcode")}
      else{
        tmp <- query %>% rownames_to_column("cell_barcode")}
      if(nrow(tmp) == 0) {
        return(data.frame())}
      cdr3 <- reference[i, ]$junction_aa_VDJ
      tmp$dist <- as.numeric(lapply(tmp$junction_aa_VDJ, function(x){stringdist::stringdist(cdr3, x, method = 'lv') / max(nchar(cdr3), nchar(x))}))
      tmp$refseq <- cdr3
      return(tmp)}}

  if(nrow(output) > 0){
    output <- output %>% 
      group_by(cell_barcode) %>%
      slice_min(n = 1, order_by = dist, with_ties = F) %>%
      ungroup() %>%
      mutate(sim = ifelse(dist <= 1, 1-dist, 0)) %>%
      dplyr::select(-c("dist")) %>%
      arrange(desc(sim))
    return(output)}
  else{
    stop("No matching sequences")
  }
}
```