## detectBCR
An R package to detect and visualize antigen-specific BCRs by cross-referencing individual query sequence to our manually curated collection of known antigen-binding BCR contigs.

The package consist 2 main functions :  

`get_reference()` allows users to request a reference database of BCR contigs specific for the antigen of interest.

`find_publicBCR()` allows users to
* Select query BCR contigs with matching heavy/light chain VJ gene usage with reference BCR contigs
* Determine how similar each query CDR3 amino acid (AA) sequence is to all antigen-specific CDR3 AA sequences

## Installation
```{r}
install.packages("devtools")
devtools::install_github("hungms/detectBCR", dependencies = T)
```


## Running the pipeline

In this tutorial we will demonstrate how we can determine **PUBLIC** Sar-CoV-2 specific BCR sequences from the example data in the package.

### Example query data
Here we will load the example QUERY data included in the package, which contain single-cell BCR sequences from plasma cells after SARS-CoV-2 mRNA-1273 vaccine. 

Each row of the data should represent a cell / unique pair of BCR sequence. The data **must** contain a column for heavy chain CDR3 AA sequence, and optionally other columns containing information about V gene, J gene, CDR3 AA sequence for heavy/light chain.

```{r}
# get example query data
path <- system.file("extdata", "dandelion_metadata.csv", package = "detectBCR")
query <- read.csv(path, row.names = 1)[1:10,] # reduce data for run time
head(query)
```

### Request reference data

```{r}
# get reference Sars-CoV2 BCR database
reference <- get_reference(antigen = "Sars-CoV-2",  org = "human")
head(reference)
```


### Running find_publicBCR()
Below we will match the query BCR contigs with the Sars-CoV reference BCR contigs



Distance calculation
```{r}
# run find_publicBCR
output <- match_CDR3(
  query = query,                                 # query dataframe
  reference = reference,                         # reference dataframe
  dist_method = c("hamming", "levenshtein"),     # CDR3 distance calculation method
  heavyCDR3 = "junction_aa_VDJ",                 # Required; query column name that stores heavyCDR3 amino acid sequence
  heavyV = "v_vall_VDJ_main",                    # Optional; query column name that stores heavyCDR3
  heavyJ = "j_call_VDJ_main",                    # Optional; query column name that stores heavyCDR3
  # lightCDR3 = NA,                                # Optional; query column name that stores lightCDR3 amino acid sequence
  # lightV = NA,                                   # Optional; query column name that stores light chain V gene usage
  # lightJ = NA,                                   # Optional; query column name that stores light chain J gene usage
  # output_dir = NULL,                             # Optional; directory to store output file
  # ncores = 1                                     # Optional; number of cores for parallel processing
)
```

### Output & downstream analysis

Running the pipeline should return a dataframe of query contigs matching with the closest antigen-specifc contigs. The dataframe should contain columns as follows:

`ref_*` = reference BCR metadata.
`dist_method` = method to calculate CDR3 amino acid distance, either `hamming` or `levenshtein`.
`dist` = CDR3 amino acid distance calculated from `dist_method`, containing values ranging from 0 to 1. 0 means the CDR3 between query and reference sequence is an exact match.
`mean_dist` = mean CDR3 amino acid distance containing values ranging from 0 to 1. This is an average of `dist` when both `heavyCDR3` and `lightCDR3` contigs are matched.

```{r}
head(output)
```


Users can define "antigen-specificity" base on their choice of method and threshold. For example, we have defined antigen-speficic BCRs base on levenshtein distance less than 0.2.

```{r}
ag_specific_bcr <- output %>%
  filter(dist_method == "levenshtein") %>%
  filter(dist < 0.2)

head(ag_specific_bcr)
dim(ag_specific_bcr)
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

## Citation
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