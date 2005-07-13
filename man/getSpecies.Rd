\name{getSpecies}
\alias{getSpecies}
\title{Displays species identifiers currently present in Ensembl}
\description{Displays species identifiers currently present in Ensembl}
\usage{getSpecies(mart, db = c("ensembl"))}

\arguments{
\item{db}{database to query, only works for Ensembl}
\item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
 

}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

if(interactive()){
mart <- martConnect()

species =getSpecies(mart = mart)
show(species)

martDisconnect(mart = mart)
}
}
\keyword{methods}

