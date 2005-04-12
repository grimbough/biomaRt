\name{getGene}
\alias{getGene}
\title{Retrieves Gene information}
\description{This function function retrieves gene information from Ensemble given a gene identifier}
\usage{getGene( id = NULL, type = NULL, array = NULL, species = NULL, db = "ensembl", mart = NULL, output = "martTable" )}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy, entrezgene, embl or refseq}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{species}{species, only to be used if type equals locuslink.  Accepted values for species can be retrieved with the function getSpecies()} 
  \item{db}{database you want to query either ensembl (default) or vega}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
  \item{output}{desired output type, either object of class martTable or GeneObject}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

mart <- martConnect()

#example using affy id from Ensembl

getGene( id = "1939_at", array = "hg_u95av2", mart = mart)

#example using locuslink id from vega

getGene( id = 100, type = "entrezgene", species = "hsapiens", db = "vega", mart = mart)

martDisconnect(mart = mart)
}
\keyword{methods}

