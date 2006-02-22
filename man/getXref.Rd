\name{getXref}
\alias{getXref}
\title{Get cross reference in Ensembl}
\description{This powerful function retrieves any cross reference in Ensembl}
\usage{getXref(id, from.species, to.species, from.xref, to.xref, db = "ensembl", mart)}

\arguments{
  \item{id}{any identifier}
  \item{from.xref}{type of identifier, value should correspond to results of getPossibleXrefs function}
  \item{to.xref}{type of identifier, value should correspond to results of getPossibleXrefs function}
  \item{from.species}{species of which the id is from.  Accepted values for species can be retrieved with the function getSpecies()} 
  \item{to.species}{species in which you want to find the homologs.  Accepted values for species can be retrieved with the function getSpecies()} 
  \item{db}{database to use. Currently this function only works with Ensembl}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Sean Davis, Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

if(interactive()){
mart <- useMart("ensembl",mysql=TRUE)

#starting from an affy id of chip hg_u95av2 and id 1939_at, give me corresponding affy identifiers on the affy mouse chip mouse430_2


xref = getXref(id="1939_at", from.species="hsapiens", to.species = "mmusculus", from.xref ="affy_hg_u95av2", to.xref="affy_mouse430_2",mart=mart)
show(xref)

martDisconnect(mart = mart)
}
}
\keyword{methods}

