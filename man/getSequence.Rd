\name{getSequence}
\alias{getSequence}
\title{Retrieves sequences}
\description{This function retrieves sequences given the chomosome, start and end position}
\usage{getSequence( species = NULL, chromosome = NULL, start = NULL, end = NULL, martTable = NULL, mart = NULL )}

\arguments{
  \item{species}{Species information.  A list of possible values for this argument can be obtained by the function getSpecies.}
  \item{chromosome}{Chromosome name}
  \item{start}{start position of sequence on chromosome}
  \item{end}{end position of sequence on chromosome}
  \item{martTable}{alternatively a martTable can be used containing all information}
  \item{mart}{object of class Mart, representing a connection to a BioMart database}
  

}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{


mart <- martConnect(mart = "sequence")

getSequence(species="gallus_gallus", chromosome = 1, start = 400, end = 500, mart = mart)

martDisconnect(mart)

}
\keyword{methods}

