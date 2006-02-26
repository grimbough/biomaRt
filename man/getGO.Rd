\name{getGO}
\alias{getGO}
\title{Retrieves GO information}
\description{This function retrieves GO information from Ensemble given a gene identifier}
\usage{getGO( id, type, array,  mart)}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy, entrezgene, embl, hugo, ensembl, unigene, agilentprobe or refseq}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\examples{
if(interactive()){

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#example using affy id

go = getGO( id = "1939_at", array = "affy_hg_u95av2", mart = mart)
show(go)

#example using locuslink id

go = getGO( id = 672, type = "entrezgene", mart = mart)
show(go)

martDisconnect(mart = mart)
}
}
\keyword{methods}

