\name{getOMIM}
\alias{getOMIM}
\title{Retrieves OMIM information}
\description{This function retrieves OMIM information from Ensemble given a gene identifier}
\usage{getOMIM( id = NULL, type = NULL, array = NULL, mart = NULL)}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy or locuslink}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{mart}{object of class Mart, representing a connection to a BioMart database}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

mart <- martConnect("ensembl")

#example using affy id

getOMIM( id = "1939_at", array = "hgu95av2", mart = mart)

#example using locuslink id

getOMIM( id = 672, type = "locuslink", mart = mart)

martDisconnect(mart)
}
\keyword{methods}

