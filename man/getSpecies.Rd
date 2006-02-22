\name{getSpecies}
\alias{getSpecies}
\title{Displays species identifiers currently present in Ensembl}
\description{Displays species identifiers currently present in Ensembl}
\usage{getSpecies(mart)}

\arguments{
\item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

if(interactive()){
mart <- useMart("ensembl",mysql=TRUE)

species =getSpecies(mart = mart)
show(species)

martDisconnect(mart = mart)
}
}
\keyword{methods}

