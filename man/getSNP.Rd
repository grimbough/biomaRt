\name{getSNP}
\alias{getSNP}
\title{Retrieves SNP's between a given chromosome start and end position}
\description{Retrieves SNP's between a given chromosome start and end position}
\usage{getSNP( species = NULL, chromosome = NULL, start = NULL, end = NULL, mart = NULL)}

\arguments{
  \item{species}{Species information.  A list of possible values for this argument can be obtained by the function getSpecies.}
  \item{chromosome}{Chromosome name}
  \item{start}{start position on given chromosome}
  \item{end}{end position on given chromosome}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

mart <- martConnect()

getSNP(chromosome = 8, start = 148350, end = 148612, species ="homo_sapiens", mart = mart)

martDisconnect(mart = mart)
}
\keyword{methods}

