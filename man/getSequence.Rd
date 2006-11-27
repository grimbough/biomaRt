\name{getSequence}
\alias{getSequence}
\title{Retrieves sequences}
\description{This function retrieves sequences given the chomosome, start and end position.  For genomic sequences using biomaRt in MySQL mode the sequence that is returned is the 5' to 3'+ strand. So if the sequence of interest is the minus strand, one has to compute the reverse complement of the retrieved sequence, which can be done using functions provided in the matchprobes package.  Using getSequence in web service mode (default) generates 5' to 3' sequences of the requested type on the correct strand.}
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

}
}
\keyword{methods}

