\name{useDataset}
\alias{useDataset}
\title{Select a dataset to use and updates Mart object}
\description{This function selects a dataset and updates the Mart object}
\usage{useDataset(dataset,mart)}

\arguments{
\item{dataset}{Dataset you want to use.  List of possible datasets can be retrieved using the function listDatasets}
\item{mart}{Mart object created with the useMart function}
}

\author{Steffen Durinck, http://www.stat.berkeley.edu/~steffen} 

\examples{
if(interactive()){
mart=useMart("ensembl")
mart=useDataset("hsapiens_gene_ensembl", mart = mart)
}
}
\keyword{methods}

