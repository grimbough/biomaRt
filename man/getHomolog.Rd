\name{getHomolog}
\alias{getHomolog}
\title{Get homologs}
\description{This function retrieves homologs of genes from one species in another species}
\usage{getHomolog(id = NULL, fromtype = NULL, totype = NULL, fromspecies = NULL, tospecies = NULL, mart = NULL, output = "martTable")}

\arguments{
  \item{id}{gene identifier}
  \item{fromtype}{type of identifier, value should be either entrezgene, embl or refseq}
  \item{totype}{type of identifier, value should be either entrezgene, embl or refseq}
  \item{fromspecies}{species of which the id is from.  Accepted values for species can be retrieved with the function getSpecies()} 
  \item{tospecies}{species in which you want to find the homologs.  Accepted values for species can be retrieved with the function getSpecies()} 
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
  \item{output}{desired output type, object of class martTable}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

mart <- martConnect()

getHomolog(id = 1:20, fromspecies = 'homo_sapiens', tospecies = 'mus_musculus', fromtype = 'entrezgene', totype = 'refseq', mart = mart)

martDisconnect(mart = mart)
}
\keyword{methods}

