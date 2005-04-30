\name{getHomolog}
\alias{getHomolog}
\title{Get homologs}
\description{This function retrieves homologs of genes from one species in another species}
\usage{getHomolog(id = NULL, from.type = NULL, to.type = NULL, from.species = NULL, to.species = NULL, mart = NULL, output = "martTable")}

\arguments{
  \item{id}{gene identifier}
  \item{from.type}{type of identifier, value should be either entrezgene, embl or refseq}
  \item{to.type}{type of identifier, value should be either entrezgene, embl or refseq}
  \item{from.species}{species of which the id is from.  Accepted values for species can be retrieved with the function getSpecies()} 
  \item{to.species}{species in which you want to find the homologs.  Accepted values for species can be retrieved with the function getSpecies()} 
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
  \item{output}{desired output type, object of class martTable}
}

\author{Sean Davis, Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
\dontrun{
mart <- martConnect()

getHomolog(id = 1:20, from.species = 'hsapiens', to.species = 'mmusculus', from.type = 'entrezgene', to.type = 'refseq', mart = mart)

martDisconnect(mart = mart)
}
}
\keyword{methods}

