\name{getINTERPRO}
\alias{getINTERPRO}
\title{Retrieves INTERPRO protein domain information}
\description{This function retrieves INTERPRO protein domain information using Ensembl and  UNIPROT (http://www.ebi.ac.uk/uniprot) given a gene identifier}
\usage{getINTERPRO( id, type, array, species, mart)}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy, entrezgene, embl, hugo, ensembl or refseq}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{species}{species, only to be used if type equals locuslink.  Accepted values for species can be retrieved by the function getSpecies()}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){
mart <- martConnect()

#example using affy id

interpro = getINTERPRO( id = "1939_at", array = "hgu95av2", mart = mart)
show(interpro)

#example using locuslink id

interpro = getINTERPRO( id = 672, type = "entrezgene", species = "hsapiens", mart = mart)
show(interpro)

martDisconnect(mart = mart)
}
}
\keyword{methods}

