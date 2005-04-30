\name{getFeature}
\alias{getFeature}
\title{Retrieves Feature starting from a gene symbol and array identifier}
\description{This function function retrieves features representing a certain gene given a gene symbol array identifier}
\usage{getFeature( symbol = NULL, array = NULL, type = NULL, mart = NULL)}

\arguments{
  \item{symbol}{gene identifier}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{type}{type of input identifier, this is either locuslink, embl or refseq.}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}

}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
\dontrun{
mart <- martConnect()

#example using affy id from Ensembl

getFeature( symbol = "P53", array = "hg_u95av2", mart = mart)

martDisconnect(mart = mart)
}
}
\keyword{methods}

