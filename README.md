## detectBCR
An R package to detect and visualize antigen-specific BCRs and clonotypes by cross-referencing individual query sequence to our manually curated collection of known antigen-binding BCR sequences.

## Functionality
* Match individual query BCR to antigen-binding BCR base on heavy/light chain VJ gene usage
* Determine minimum `levenshtein` or `hamming` distances for query CDR3 sequences to our antigen-binding BCR sequences
* Determine convergent, antigen-specific clones

## Installation
```{r}
install.packages("devtools")
devtools::install_github("hungms/detectBCR", dependencies = T)
```

## Example Run
```{r}
# load example data
path <- system.file("extdata", "dandelion_metadata.csv", package = "detectBCR")

# reduce data for run time
query <- read.csv(path, row.names = 1)[1:100,]

# run match_CDR3
output <- match_CDR3(
  query = query,                                 # query dataframe
  antigen = "Sars-CoV-2",                        # antigen to search against - "Sars-CoV-2", "Vaccinia", "Tetanus", "Measles", "Mumps"
  org = "human",                                 # organism
  dist_method = c("hamming", "levenshtein"),     # CDR3 distance calculation method
  heavyCDR3 = "junction_aa_VDJ",                 # Required; column name that stores heavyCDR3 amino acid sequence
  heavyV = NA,                                   # Optional; column name that stores heavyCDR3
  heavyJ = NA,                                   # Optional; column name that stores heavyCDR3
  lightCDR3 = NA,                                # Optional; column name that stores lightCDR3 amino acid sequence
  lightV = NA,                                   # Optional; column name that stores light chain V gene usage
  lightJ = NA,                                   # Optional; column name that stores light chain J gene usage
  output_dir = NULL,                             # Optional; directory to store output file
  ncores = 1                                     # Optional; number of cores for parallel processing
)
```


## Collection
Currently our collection contains BCR sequences bind to the following antigens : 
* `Sars-CoV-2` 
* `Vaccinia`
* `Tetanus` - TBC
* `Measles` - TBC
* `Mumps` - TBC

```
| VERSION | ANTIGEN    | DATABASE           | REANALYSIS |  NOTE | DOI  |
| ------- | ---------- | ------------------ | ---------- | ---- | ---- |
| v0.0.0  | Sars-CoV-2 | CoV-AbDab          | N          |      | 10.1093/bioinformatics/btaa739 |
| v0.0.0  | Vaccinia   | Chappert_2022      | Y          | B5+  | 10.1016/j.immuni.2022.08.019   |
| v0.0.1  | Sars-CoV-2 | LopezDeAssis_2023  | Y          | S2P+ | 10.1016/j.celrep.2023.112780   |
| v0.0.1  | Sars-CoV-2 | FerreiraGomez_2024 | Y          | TBC  | 10.1038/s41467-024-48570-0     |
| v0.0.1  | Tetanus    | FerreiraGomez_2024 | Y          | TBC  | 10.1038/s41467-024-48570-0     |    
| v0.0.1  | Measles    |                    |            |      | 


```

## Glossary
Terms :
`Public Clones` = Clones containing BCR sequence matching our antigen-binding BCR database
`Convergent Clones` = Antigen-specific clones that are present in multiple individuals

Abbreviations :
`heavyV` = Heavy chain V gene
`heavyJ` = Heavy chain J gene
`heavyCDR3` = Heavy chain CDR3 amino acid sequence
`lightV` = Light chain V gene
`lightJ` = Light chain J gene
`ligjtCDR3` = Light chain CDR3 amino acid sequence

## Inspiration
Codes were adapted from previous publication in Cell Reports :

```
@article{LOPESDEASSIS2023112780,
title = {Tracking B cell responses to the SARS-CoV-2 mRNA-1273 vaccine},
journal = {Cell Reports},
volume = {42},
number = {7},
pages = {112780},
year = {2023},
issn = {2211-1247},
doi = {https://doi.org/10.1016/j.celrep.2023.112780},
url = {https://www.sciencedirect.com/science/article/pii/S221112472300791X},
author = {Felipe {Lopes de Assis} and Kenneth B. Hoehn and Xiaozhen Zhang and Lela Kardava and Connor D. Smith and Omar {El Merhebi} and Clarisa M. Buckner and Krittin Trihemasava and Wei Wang and Catherine A. Seamon and Vicky Chen and Paul Schaughency and Foo Cheung and Andrew J. Martins and Chi-I Chiang and Yuxing Li and John S. Tsang and Tae-Wook Chun and Steven H. Kleinstein and Susan Moir},
keywords = {SARS-CoV-2, mRNA vaccine, B cells, immunological memory, single-cell profiling, BCR repertoire},
abstract = {Summary
Protective immunity following vaccination is sustained by long-lived antibody-secreting cells and resting memory B cells (MBCs). Responses to two-dose SARS-CoV-2 mRNA-1273 vaccination are evaluated longitudinally by multimodal single-cell analysis in three infection-naïve individuals. Integrated surface protein, transcriptomics, and B cell receptor (BCR) repertoire analysis of sorted plasmablasts and spike+ (S-2P+) and S-2P− B cells reveal clonal expansion and accumulating mutations among S-2P+ cells. These cells are enriched in a cluster of immunoglobulin G-expressing MBCs and evolve along a bifurcated trajectory rooted in CXCR3+ MBCs. One branch leads to CD11c+ atypical MBCs while the other develops from CD71+ activated precursors to resting MBCs, the dominant population at month 6. Among 12 evolving S-2P+ clones, several are populated with plasmablasts at early timepoints as well as CD71+ activated and resting MBCs at later timepoints, and display intra- and/or inter-cohort BCR convergence. These relationships suggest a coordinated and predictable evolution of SARS-CoV-2 vaccine-generated MBCs.}
}
```