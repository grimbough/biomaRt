\name{getGO}
\alias{getGO}
\title{Retrieves GO information}
\description{This function retrieves GO information from Ensemble given a gene identifier}
\usage{getGO( id = NULL, type = NULL, array = NULL, species = NULL, mart = NULL)}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy or locuslink}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{species}{species, only to be used if type equals locuslink.  Accepted values for species can be retrieved by the function getSpecies()}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

mart <- martConnect()

#example using affy id

getGO( id = "1939_at", array = "hgu95av2", mart = mart)

#example using locuslink id

getGO( id = 672, type = "locuslink", species = "homo_sapiens", mart = mart)

martDisconnect(mart = mart)
}
\keyword{methods}

