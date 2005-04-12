\name{getSequence}
\alias{getSequence}
\title{Retrieves sequences}
\description{This function retrieves sequences given the chomosome, start and end position}
\usage{getSequence( species = NULL, chromosome = NULL, start = NULL, end = NULL, martTable = NULL, mart = NULL)}

\arguments{
  \item{species}{Species information.  A list of possible values for this argument can be obtained by the function getSpecies.}
  \item{chromosome}{Chromosome name}
  \item{start}{start position of sequence on chromosome}
  \item{end}{end position of sequence on chromosome}
  \item{martTable}{alternatively a martTable can be used containing all information}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

mart <- martConnect()

getSequence(species="ggallus", chromosome = 1, start = 400, end = 500, mart = mart)

martDisconnect(mart = mart)

}
\keyword{methods}

