\name{getFeature}
\alias{getFeature}
\title{Retrieves Feature starting from a gene symbol and array identifier}
\description{This function function retrieves features representing a certain gene given a gene symbol array identifier}
\usage{getFeature( symbol = NULL, array = NULL, type = NULL, mart = NULL )}

\arguments{
  \item{symbol}{gene identifier}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{type}{type of input identifier, this is either locuslink, embl or refseq.}
  \item{mart}{object of class Mart, representing a connection to a BioMart database}
  

}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

#example using affy id from Ensembl

mart <- martConnect(mart = "ensembl")

getFeature( symbol = "P53", array = "hgu95av2", mart = mart)

martDisconnect(mart)

}
\keyword{methods}

