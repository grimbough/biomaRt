\name{getGO}
\alias{getGO}
\title{Retrieves GO information}
\description{This function retrieves GO information from Ensembl given a gene identifier}
\usage{getGO( id, type, mart)}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, possible values can be obtained by the listFilters function.  Examples are  entrezgene, hgnc\_symbol (for hugo gene sy
mbol), ensembl\_gene\_id, unigene, agilentprobe, affy\_hg\_u133\_plus\_2,  refseq\_dna, etc.}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#example using affy id

go = getGO( id = "1939_at", type = "affy_hg_u95av2", mart = mart)
show(go)

#example using entrezgene id

go = getGO( id = 672, type = "entrezgene", mart = mart)
show(go)
}
}
\keyword{methods}

