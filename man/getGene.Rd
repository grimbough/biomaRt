\name{getGene}
\alias{getGene}
\title{Retrieves gene annotation information}
\description{This function function retrieves gene annotations (chromsome name, band, start position, end position, gene description and gene symbol) from Ensembl given many types of  identifiers}
\usage{getGene( id, type, mart)}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, possible values can be obtained by the listFilters function.  Examples are  entrezgene, hgnc\_symbol (for hugo gene symbol), ensembl\_gene\_id, unigene, agilentprobe, affy\_hg\_u133\_plus\_2,  refseq\_dna, etc.}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

#example using affy id

g = getGene( id = "1939_at", type = "affy_hg_u95av2", mart = mart)
show(g)

#example using Entrez Gene id

g = getGene( id = "100", type = "entrezgene", mart = mart)
show(g)
}
}
\keyword{methods}

