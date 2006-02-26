\name{getHomolog}
\alias{getHomolog}
\title{Get homologs}
\description{This function retrieves homologs of genes from one species in another species}
\usage{getHomolog(id, from.type, to.type, from.array, to.array, from.mart, to.mart)}

\arguments{
  \item{id}{gene identifier}
  \item{from.type}{type of identifier, value should be either entrezgene, embl or refseq}
  \item{to.type}{type of identifier, value should be either entrezgene, embl or refseq}
  \item{from.array}{Affy array from which the input identifiers come.  Note if this has a value then you don't have to specify from.type and from.species.}
  \item{to.array}{Affy array to which you want to map.  Note if this has a value then you don't have to specify to.type and to.species.}
  \item{from.mart}{Mart object using dataset from species of which the query id is from.} 
  \item{to.mart}{Mart object using dataset from species you want to find the homologs off.} 
}

\author{Sean Davis, Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\examples{
if(interactive()){
from.mart <- useMart("ensembl","hsapiens_gene_ensembl")
to.mart <- useMart("ensembl","mmusculus_gene_ensembl")

#HUGO to Entrez Gene

homolog = getHomolog(id = 1:20, from.mart = from.mart, to.mart = to.mart, from.type = 'entrezgene', to.type = 'refseq')
show(homolog)

#ensembl to ensembl

homolog = getHomolog( id = "ENSG00000072778", from.mart = from.mart, from.type = "ensembl", to.type="ensembl",to.mart=to.mart)
show(homolog)

#Affy to Affy

homolog = getHomolog( id = "1424184_at", to.array = "hg_u95av2", from.array = "mouse430a_2", from.mart = from.mart, to.mart=to.mart )
show(homolog)

#Affy to ensembl

homolog = getHomolog( id = "1424184_at", to.type = "ensembl", from.array = "mouse430a_2", from.mart = from.mart, to.mart = to.mart )
show(homolog)

}
}
}
\keyword{methods}

