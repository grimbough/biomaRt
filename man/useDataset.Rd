\name{useDataset}
\alias{useDataset}
\title{Select a dataset to use and updates Mart object}
\description{This function selects a dataset and updates the Mart object}
\usage{useDataset(dataset,mart, verbose = FALSE)}

\arguments{
\item{dataset}{Dataset you want to use.  List of possible datasets can be retrieved using the function listDatasets}
\item{mart}{Mart object created with the useMart function}
\item{verbose}{Give detailed output of what the method is doing, for debugging}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
mart=useMart("ensembl")
mart=useDataset("hsapiens_gene_ensembl", mart = mart)
}
}
\keyword{methods}

