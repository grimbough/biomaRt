\name{listFilters}
\alias{listFilters}
\title{lists the filters available in the selected dataset}
\description{lists the filters available in the selected dataset}
\usage{listFilters(mart)}

\arguments{
\item{mart}{object of class Mart, containing connections to the BioMart database.  You have to create such an object using the function useMart and update this Mart object by selecting a dataset.}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){
#marts <- listMarts()
#index<-grep("ensembl",marts)

#mart <- useMart(marts[index])
#datasets <- listDatasets(mart)
#mart <- useDataset(datasets[1],mart = mart)
#filters <- listFilters(mart)
#filters[1:10]
#martDisconnect(mart = mart)
}
}
\keyword{methods}

