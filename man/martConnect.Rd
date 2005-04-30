\name{martConnect}
\alias{martConnect}
\title{Connects to available BioMarts}
\description{This function connects to available BioMarts}
\usage{martConnect(host = "ensembldb.ensembl.org", user = "anonymous", password = "")}

\arguments{
\item{host}{host where database resides}
\item{user}{username}
\item{password}{password}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
\dontrun{
mart <- martConnect()

martDisconnect(mart = mart)
}
}
\keyword{methods}

