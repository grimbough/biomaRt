\name{getFeature}
\alias{getFeature}
\title{Retrieves Feature starting from a gene symbol and array identifier}
\description{This function function retrieves features representing a certain gene given a gene symbol array identifier}
\usage{getFeature( symbol, OMIM, OMIMID, GO, GOID, array, species, chromosome, start, end, type, mart)}

\arguments{
  \item{symbol}{gene symbol}
  \item{OMIM}{OMIM term e.g. diabetes} 
  \item{OMIMID}{OMIM identifier} 
  \item{GO}{GO term}
  \item{GOID}{GO identifier}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{species}{species}
  \item{chromosome}{chromosome, can be used without start and end positions}
  \item{start}{start on chromosome}
  \item{end}{end on chromosome}
  \item{type}{type of input identifier, this is either locuslink, embl or refseq.}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}

}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){
mart <- useMart("ensembl",mysql=TRUE)

#example using affy id from Ensembl

feat = getFeature( symbol = "P53", array = "hg_u95av2", mart = mart)
show(feat)

martDisconnect(mart = mart)
}
}
\keyword{methods}

