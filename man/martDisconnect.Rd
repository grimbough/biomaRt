\name{martDisconnect}
\alias{martDisconnect}
\title{Disconnects from BioMarts}
\description{This function closes connections of the Mart object the BioMarts}
\usage{martConnect()}

\arguments{
\item{mart}{Mart object created with function martConnect}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

mart <- martConnect()

martDisconnect(mart = mart)
}
\keyword{methods}

