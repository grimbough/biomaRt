\name{martDisconnect}
\alias{martDisconnect}
\title{Disconnects from BioMarts}
\description{This function closes connections of the Mart object the BioMarts}
\usage{martDisconnect(mart)}

\arguments{
\item{mart}{Mart object created with function martConnect}
}

\author{Steffen Durinck, \url{http://www.stat.berkeley.edu/~steffen}} 

\examples{
if(interactive()){
mart <- martConnect()

martDisconnect(mart = mart)
}
}
\keyword{methods}

