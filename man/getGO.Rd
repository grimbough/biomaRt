\name{getGO}
\alias{getGO}
\title{Retrieves GO information}
\description{This function retrieves GO information from Ensemble given a gene identifier}
\usage{getGO( id, type, array, species, mart, output="data.frame")}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy, entrezgene, embl, hugo, ensembl or refseq}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{species}{species, only to be used if type equals locuslink.  Accepted values for species can be retrieved by the function getSpecies()}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
  \item{output}{sepcify output format, default is data.frame.  You can also select a martTable as output}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){
mart <- martConnect()

#example using affy id

go = getGO( id = "1939_at", array = "hg_u95av2", mart = mart)
show(go)

#example using locuslink id

go = getGO( id = 672, type = "entrezgene", species = "hsapiens", mart = mart)
show(go)

martDisconnect(mart = mart)
}
}
\keyword{methods}

