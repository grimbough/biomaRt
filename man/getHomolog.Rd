\name{getHomolog}
\alias{getHomolog}
\title{Get homologs}
\description{This function retrieves homologs of genes from one species in another species}
\usage{getHomolog(id, from.type, to.type, from.mart, to.mart)}

\arguments{
  \item{id}{gene identifier}
  \item{from.type}{type of identifier of the input. Possible values depend on the selected species and can be obtained by the listFilters function and using the from.mart.  Examples for hsapiens are: entrezgene, hgnc\_symbol (for hugo gene symbol), ensembl\_gene\_id, unigene, agilentprobe, affy\_hg\_u133\_plus\_2,  refseq\_dna, etc.}
  \item{to.type}{type of identifier that needs to be retrieved from as homolog id. Possible values depend on the selected species and can be obtained by the listAttributes function and using the to.mart.  Examples for hsapiens are entrezgene, hgnc\_symbol (for hugo gene symbol), ensembl\_gene\_id, unigene, agilentprobe, affy\_hg\_u133\_plus\_2,  refseq\_dna, etc. }
  \item{from.mart}{Mart object using dataset from species of which the query id is from.} 
  \item{to.mart}{Mart object using dataset from species you want to find the homologs off.} 
}

\author{Steffen Durinck, Sean Davis} 

\examples{
if(interactive()){
from.mart <- useMart("ensembl","hsapiens_gene_ensembl")
to.mart <- useMart("ensembl","mmusculus_gene_ensembl")

#HUGO to Entrez Gene

homolog = getHomolog(id = 1:20, from.mart = from.mart, to.mart = to.mart, from.type = 'entrezgene', to.type = 'refseq_dna')
show(homolog)

#ensembl to ensembl

homolog = getHomolog( id = "ENSG00000072778", from.mart = from.mart, from.type = "ensembl_gene_id", to.type="ensembl_gene_id",to.mart=to.mart)
show(homolog)

#Affy to Affy

homolog = getHomolog( id = "1939_at", to.type = "affy_mouse430_2", from.type = "affy_hg_u95av2", from.mart = from.mart, to.mart=to.mart )
show(homolog)

#Ensembl to Affy

homolog = getHomolog( id = "ENSG00000072778", to.type = "affy_mouse430_2", from.type = "ensembl_gene_id", from.mart = from.mart, to.mart = to.mart )
show(homolog)

}
}
\keyword{methods}

