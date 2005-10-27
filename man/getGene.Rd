\name{getGene}
\alias{getGene}
\title{Retrieves Gene information}
\description{This function function retrieves gene information from Ensemble given a gene identifier}
\usage{getGene( id, type, array, species, db = "ensembl", mart, output="data.frame")}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy, entrezgene, embl, hugo, ensembl or refseq}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{species}{species, only to be used if type equals locuslink.  Accepted values for species can be retrieved with the function getSpecies()} 
  \item{db}{database you want to query either ensembl (default) or vega}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
  \item{output}{Specify output, default is a data.frame. You can also select a martTable as output but this is going to dissappear}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

if(interactive()){

mart <- martConnect( biomarts = c("ensembl","vega"))

#example using affy id from Ensembl

g = getGene( id = "1939_at", array = "hg_u95av2", mart = mart)
show(g)

#example using locuslink id from vega

g = getGene( id = 100, type = "entrezgene", species = "hsapiens", db = "vega", mart = mart)
show(g)

martDisconnect(mart = mart)

}
}
}
\keyword{methods}

