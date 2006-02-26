\name{getGene}
\alias{getGene}
\title{Retrieves Gene information}
\description{This function function retrieves gene information from Ensemble given a gene identifier}
\usage{getGene( id, type, array, mart)}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either entrezgene, embl, hugo, ensembl, unigene, agilentprobe or refseq}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\examples{

if(interactive()){

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

#example using affy id

g = getGene( id = "1939_at", array = "affy_hg_u95av2", mart = mart)
show(g)

#example using Entrez Gene id

g = getGene( id = "100", type = "entrezgene", mart = mart)
show(g)

martDisconnect(mart = mart)

}
}
}
\keyword{methods}

