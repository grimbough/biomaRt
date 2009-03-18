\name{martDisconnect}
\alias{martDisconnect}
\title{Disconnects from BioMarts (defunct)}
\description{The function 'martDisconnect' is defunct. In BioMart webservice mode, disconnecting is not necessary. Background explanation: In previous versions of this package, a MySQL mode was also supported, which required disconnecting.}

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

