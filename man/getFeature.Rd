\name{getFeature}
\alias{getFeature}
\title{Retrieves identifiers starting from a gene symbol, chromosomal position or GO identifier}
\description{This function function retrieves identifiers representing a certain gene given a gene symbol, chromosomal positions or GO identifier}
\usage{getFeature( symbol, OMIMID, GOID, chromosome, start, end, type, mart)}

\arguments{
  \item{symbol}{gene symbol}
  \item{OMIMID}{OMIM identifier} 
  \item{GOID}{GO identifier}
  \item{chromosome}{chromosome, can be used without start and end positions}
  \item{start}{start on chromosome}
  \item{end}{end on chromosome}
  \item{type}{type of identifier that needs to be retrieved, possible values can be obtained by the listFilters function.  Examples are  entrezgene, hgnc\_symbol (for hugo gene symbol), ensembl\_gene\_id, unigene, agilentprobe, refseq\_dna, affy\_hg\_u133\_plus\_2,  etc.}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function martConnect.}

}

\author{Steffen Durinck} 

\examples{
if(interactive()){
mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

#example using affy id from Ensembl

feat = getFeature( symbol = "TP53",type  = "affy_hg_u95av2", mart = mart)
show(feat)
}
}
\keyword{methods}

