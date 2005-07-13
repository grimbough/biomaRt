\name{getHomolog}
\alias{getHomolog}
\title{Get homologs}
\description{This function retrieves homologs of genes from one species in another species}
\usage{getHomolog(id, from.type, to.type, from.array, to.array, from.species, to.species, mart, output = "martTable")}

\arguments{
  \item{id}{gene identifier}
  \item{from.type}{type of identifier, value should be either entrezgene, embl or refseq}
  \item{to.type}{type of identifier, value should be either entrezgene, embl or refseq}
  \item{from.array}{Affy array from which the input identifiers come.  Note if this has a value then you don't have to specify from.type and from.species.}
  \item{to.array}{Affy array to which you want to map.  Note if this has a value then you don't have to specify to.type and to.species.}
  \item{from.species}{species of which the id is from.  Accepted values for species can be retrieved with the function getSpecies()} 
  \item{to.species}{species in which you want to find the homologs.  Accepted values for species can be retrieved with the function getSpecies()} 
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}
  \item{output}{desired output type, object of class martTable}
}

\author{Sean Davis, Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){
mart <- martConnect()

#HUGO to Entrez Gene

homolog = getHomolog(id = 1:20, from.species = 'hsapiens', to.species = 'mmusculus', from.type = 'entrezgene', to.type = 'refseq', mart = mart)
show(homolog)

#ensembl to ensembl

homolog = getHomolog( id = "ENSG00000072778", from.species = "hsapiens", from.type = "ensembl", to.type="ensembl",to.species="mmusculus", mart = mart )
show(homolog)

#Affy to Affy

homolog = getHomolog( id = "1424184_at", to.array = "hg_u95av2", from.array = "mouse430a_2", mart = mart )
show(homolog)

#Affy to ensembl

homolog = getHomolog( id = "1424184_at", to.species = "hsapiens", to.type = "ensembl", from.array = "mouse430a_2", mart = mart )
show(homolog)

martDisconnect(mart = mart)
}
}
}
\keyword{methods}

