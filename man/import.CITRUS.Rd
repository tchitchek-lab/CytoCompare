% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R
\name{import.CITRUS}
\alias{import.CITRUS}
\title{Importation of cell cluster profiles from a Citrus result}
\usage{
import.CITRUS(file, dictionary = NULL, exclude = c("Time",
  "fileEventNumber", "Cell_length", "fileId"), bin.width = 0.05,
  minimumClusterSizePercent = 0.05, cluster.selection = NULL,
  arcsinh.transform = FALSE)
}
\arguments{
\item{file}{a character indicating the location of the citrusClustering.Rdata file}

\item{dictionary}{a two-column data.frame providing the correspondence between the original marker names (first column) and the new marker names (second column)}

\item{exclude}{a vector containing the marker names to be excluded in the import procedure}

\item{bin.width}{a numeric value indicating the width of the bins for the marker expression densities computations}

\item{minimumClusterSizePercent}{a numeric value indicating the minimal ratio of cells per cluster to import}

\item{cluster.selection}{a character vector containing the names of the clusters to import}

\item{arcsinh.transform}{a logical value indicating if Citrus expression values must be arcsinh transformed}
}
\value{
a S4 object of class CLUSTER
}
\description{
Imports one or several cell cluster profiles identified by the Citrus algorithm into a CLUSTER object.
}
\details{
Citrus is an algorithm that clusters cells using a hierarchical clustering procedure (similarly to SPADE) and then identifies the cell clusters that are significantly associated with different biological condition phenotypes (PMID:24979804).
}
