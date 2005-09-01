\name{listAttributes}
\alias{listAttributes}
\title{lists the attributes available in the selected dataset}
\description{lists the datasets available in the selected dataset}
\usage{listAttributes(mart)}

\arguments{
\item{mart}{object of class Mart, containing connections to the BioMart database.  You have to create such an object using the function useMart and update this Mart object by selecting a dataset.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

if(interactive()){

marts <- listMarts()
index<-grep("ensembl",marts)
mart <- useMart(marts[index])
datasets <- listDatasets(mart)
mart <- useDataset(datasets[1],mart = mart)
attributes <- listAttributes(mart)
attributes[1:10]
martDisconnect(mart = mart)
}
}
\keyword{methods}

