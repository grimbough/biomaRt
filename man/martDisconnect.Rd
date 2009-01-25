\name{martDisconnect}
\alias{martDisconnect}
\title{Disconnects from BioMarts}
\description{This function closes connections of the Mart object the BioMarts and only has to be used when one connects via MySQL.  As this package only supports web service based queries this function is defunct.}
\usage{martDisconnect(mart)}

\arguments{
\item{mart}{Mart object created with the useMart function}
}

\author{Steffen Durinck, \url{http://www.stat.berkeley.edu/~steffen}} 

\examples{
if(interactive()){
mart <- martConnect()

martDisconnect(mart = mart)
}
}
\keyword{methods}

