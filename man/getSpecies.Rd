\name{getSpecies}
\alias{getSpecies}
\title{Displays species identifiers supported by the Ensembl package}
\description{Displays species identifiers supported by the Ensembl package}
\usage{getSpecies( mart = NULL)}

\arguments{
\item{mart}{object of class Mart, representing a connection to a BioMart database}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
mart <- martConnect("ensembl")
getSpecies(mart)
martDisconnect(mart)

}
\keyword{methods}

