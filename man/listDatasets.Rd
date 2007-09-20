\name{listDatasets}
\alias{listDatasets}
\title{lists the datasets available in the selected BioMart database}
\description{Lists the datasets available in the selected BioMart database}
\usage{listDatasets(mart)}

\arguments{
\item{mart}{object of class Mart created with the useMart function}
}

\author{Steffen Durinck, \url{http://www.stat.berkeley.edu/~steffen}} 

\examples{

if(interactive()){

#marts <- listMarts()
#index<-grep("ensembl",marts)

#mart <- useMart(marts[index])

#listDatasets(mart = mart)

#martDisconnect(mart = mart)
}
}
\keyword{methods}

