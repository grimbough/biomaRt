\name{getOMIM}
\alias{getOMIM}
\title{Retrieves OMIM information}
\description{This function retrieves OMIM information from Ensemble given a gene identifier}
\usage{getOMIM( id = NULL, type = NULL, array = NULL, mart = NULL)}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy, entrezgene, embl, hugo, ensembl or refseq}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
\dontrun{
mart <- martConnect()

#example using affy id

getOMIM( id = "1939_at", array = "hg_u95av2", mart = mart)

#example using locuslink id

getOMIM( id = 672, type = "entrezgene", mart = mart)

martDisconnect(mart = mart)
}
}
\keyword{methods}

