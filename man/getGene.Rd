\name{getGene}
\alias{getGene}
\title{Retrieves Gene information}
\description{This function function retrieves gene information from Ensemble given a gene identifier}
\usage{getGene( id = NULL, type = NULL, array = NULL, species = NULL, mart = NULL, output = "martTable" )}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy or locuslink}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{species}{species, only to be used if type equals locuslink.  Accepted values for species currently are: can be retrieved with the function getSpecies()} 
  \item{mart}{object of class Mart, representing a connection to a BioMart database}
  \item{output}{desired output type, either object of class martTable or GeneObject}

}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
#example using affy id from Ensembl

mart <- martConnect("ensembl")

getGene( id = "1939_at", array = "hgu95av2", mart = mart)

martDisconnect(mart)

#example using locuslink id from vega

mart<-martConnect("vega")

getGene( id = 100, type = "locuslink", species = "homo_sapiens", mart = mart)

martDisconnect(mart)
}
\keyword{methods}

