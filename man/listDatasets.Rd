\name{listDatasets}
\alias{listDatasets}
\title{lists the datasets available in the selected BioMart database}
\description{Lists the datasets available in the selected BioMart database}
\usage{listDatasets(mart, verbose = FALSE)}

\arguments{
\item{mart}{object of class Mart created with the useMart function}
\item{verbose}{Give detailed output of what the method is doing, for debugging purposes}
}

\author{Steffen Durinck} 

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

