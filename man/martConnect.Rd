\name{martConnect}
\alias{martConnect}
\title{Establish a connection with a BioMart database}
\description{This function establishes a connection with a BioMart database and returns an object of class Mart}
\usage{martConnect(mart = NULL, driver = NULL, host = NULL, dbname = NULL, user = NULL, passwd = NULL )}

\arguments{
  \item{mart}{string: name of BioMart database to connect to.  Currently supported databases are:  "ensembl", "vega"}
  \item{driver}{database driver}
  \item{host}{host of database}
  \item{dbname}{database name}
  \item{user}{username}
  \item{passwd}{password}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
#connect to Ensembl
mart <- martConnect("ensembl")

martDisconnect(mart)
}
\keyword{methods}

