% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_publicBCR.R
\name{find_publicBCR}
\alias{find_publicBCR}
\title{Detect Public BCR Sequences}
\usage{
find_publicBCR(
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
  output_dir = NULL
)
}
\arguments{
\item{query}{The query data frame}

\item{reference}{The reference data frame}

\item{heavyCDR3}{The heavy chain CDR3 sequence}

\item{heavyV}{The heavy chain V gene}

\item{heavyJ}{The heavy chain J gene}

\item{lightCDR3}{The light chain CDR3 sequence}

\item{lightV}{The light chain V gene}

\item{lightJ}{The light chain J gene}

\item{dist_method}{The distance method to use}

\item{ncores}{The number of cores to use}

\item{output_dir}{The output directory}
}
\value{
A data frame with the matched CDR3 sequences
}
\description{
Detect Public BCR Sequences
}
