\name{useDataset}
\alias{useDataset}
\title{Select a dataset to use and updates Mart object}
\description{This function selects a dataset and updates the Mart object}
\usage{useDataset(dataset,mart)}

\arguments{
\item{dataset}{Dataset you want to use.  List of possible datasets can be retrieved using the function listDatasets}
\item{mart}{Mart object created with the function useMart}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\examples{
if(interactive()){

#marts <- listMarts()
#index<-grep("ensembl",marts)

#mart <- useMart(marts[index])
#datasets <- listDatasets(mart=mart)
#mart<-useDataset(datasets[1],mart=mart)
#martDisconnect(mart = mart)
}
}
\keyword{methods}

