\name{getOMIM}
\alias{getOMIM}
\title{Retrieves OMIM information}
\description{This function retrieves OMIM information from Ensemble given a gene identifier}
\usage{getOMIM( id, type, array, mart)}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy, entrezgene, embl, hugo, ensembl, unigene, agilentprobe or refseq}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#example using affy id

omim = getOMIM( id = "1939_at", array = "affy_hg_u95av2", mart = mart)
show(omim)

#example using locuslink id

omim = getOMIM( id = 672, type = "entrezgene", mart = mart)
show(omim)

martDisconnect(mart = mart)
}
}
\keyword{methods}

