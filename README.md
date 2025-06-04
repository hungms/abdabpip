## detectBCR
An R package to detect and visualize antigen-specific BCRs and clonotypes by cross-referencing individual query sequence to our manually curated collection of known antigen-binding BCR sequences.


## Inspiration
Codes developed from 

## Functionality
* Match individual query BCR to antigen-binding BCR base on heavy/light chain VJ gene usage
* Determine minimum `levenshtein` or `hamming` distances for query CDR3 sequences to our antigen-binding BCR sequences
* Determine convergent, antigen-specific clones

## Installation
```{r}
install.packages("devtools")
devtools::install_github("hungms/detectBCR", dependencies = T)
```


## Collection
Currently our collection contains BCR sequences bind to the following antigens : 
* `Sars-CoV-2`
* `Tetanus`
* `Vaccinia`
* `Measles`
* `Mumps`

```
| ANTIGEN    | DATABASE  | REANALYSIS | NOTE | DOI  |
| ---------- | --------- | ---------- | ---- | ---- |
| Sars-Cov-2 | CoV-AbDab | N          |      | 10.1093/bioinformatics/btaa739 |
```

## Glossary
Terms :
`Public Clones` = Clones containing BCR sequence matching our antigen-binding BCR database
`Convergent Clones` = Antigen-specific clones that are present in multiple individuals

Abbreviations :
`VDJ_V` = Heavy chain V gene
`VDJ_J` = Heavy chain J gene
`VDJ_CDR3` = Heavy chain CDR3 amino acid sequence
`VJ_V` = Light chain V gene
`VJ_V` = Light chain J gene
`VJ_CDR3` = Light chain CDR3 amino acid sequence

