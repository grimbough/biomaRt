\name{getSequence}
\alias{getSequence}
\title{Retrieves sequences}
\description{This function retrieves sequences given the chomosome, start and end position.  For genomic sequences the strand that is returned is the 5' t0 3' strand. Note that using this function with MySQL is limited to genomic sequences.}
\usage{getSequence( chromosome, start, end, id, type, seqType, mart)}

\arguments{
  \item{chromosome}{Chromosome name}
  \item{start}{start position of sequence on chromosome}
  \item{end}{end position of sequence on chromosome}
  \item{id}{ An identifier or vector of identifiers.}
  \item{type}{The type of identifier used.  Supported types are hugo, ensembl, embl, entrezgene, refseq, ensemblTrans and unigene.}
  \item{seqType}{Type of sequence that you want to retrieve.  Allowed seqTypes are: cdna, peptide, 3utr, 5utr, genomic}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\examples{
if(interactive()){
mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

seq = getSequence(id="BRCA1", type="hugo", seqType="peptide", mart = mart)
show(seq)

martDisconnect(mart = mart)
}
}
\keyword{methods}

