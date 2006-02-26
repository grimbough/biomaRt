\name{getINTERPRO}
\alias{getINTERPRO}
\title{Retrieves INTERPRO protein domain information}
\description{This function retrieves INTERPRO protein domain information using Ensembl and  UNIPROT (http://www.ebi.ac.uk/uniprot) given a gene identifier}
\usage{getINTERPRO( id, type, array, mart)}

\arguments{
  \item{id}{gene identifier}
  \item{type}{type of identifier, value should be either affy, entrezgene, embl, hugo, ensembl or refseq}
  \item{array}{array identifier only to be used if gene identifiers are from affy arrays. An overview of the accepted identifiers can be obtained with the function: getAffyArrays()}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\examples{
if(interactive()){

mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#example using affy id

interpro = getINTERPRO( id = "1939_at", array = "affy_hg_u95av2", mart = mart)
show(interpro)

#example using locuslink id

interpro = getINTERPRO( id = 672, type = "entrezgene", mart = mart)
show(interpro)

}
}
\keyword{methods}

