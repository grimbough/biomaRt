\name{listDatasets}
\alias{listDatasets}
\title{lists the datasets available in the selected BioMart database}
\description{lists the datasets available in the selected BioMart database}
\usage{listDatasets(mart)}

\arguments{
\item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
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

