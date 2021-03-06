% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{fpkm}
\alias{fpkm}
\title{Calculate FPKM/RPKM}
\usage{
fpkm(mat, geneLength, libSize=NULL)
}
\arguments{
\item{mat}{data frame count matrix with rows as genes and columns as samples}

\item{geneLength}{data frame with two columns: GeneId corresponding to rownames(mat) and geneLength with the exon lengths of each gene}

\item{libSize}{numeric vector of length ncol(mat) containing (in order) the library sizes for each sample in mat. Leave NULL to use colSums(mat).}
}
\value{
data frame with FPKM/RPKM normalized gene counts
}
\description{
Function to calculate fragments/reads per kilobase per million using a gene count matrix
}
