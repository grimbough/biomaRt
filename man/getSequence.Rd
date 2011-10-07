\name{getSequence}
\alias{getSequence}
\title{Retrieves sequences}
\description{This function retrieves sequences given the chomosome, start and end position or a list of identifiers. Using getSequence in web service mode (default) generates 5' to 3' sequences of the requested type on the correct strand. The type of sequence returned can be specified by the seqType argument which takes the following values: 'cdna';'peptide' for protein sequences;'3utr' for 3' UTR sequences,'5utr' for 5' UTR sequences; 'gene\_exon' for exon sequences only; 'transcript\_exon\_intron' gives the full unspliced transcript, that is exons + introns;'gene\_exon\_intron' gives the exons + introns of a gene;'coding' gives the coding sequence only;'coding\_transcript\_flank' gives the flanking region of the transcript including the UTRs, this must be accompanied with a given value for the upstream or downstream attribute;'coding\_gene\_flank' gives the flanking region of the gene including the UTRs, this must be accompanied with a given value for the upstream or downstream attribute; 'transcript\_flank' gives the flanking region of the transcript exculding the UTRs, this must be accompanied with a given value for the upstream or downstream attribute; 'gene\_flank' gives the flanking region of the gene excluding the UTRs, this must be accompanied with a given value for the upstream or downstream attribute.
In MySQL mode the getSequence function is more limited and the sequence that is returned is the 5' to 3'+ strand of the genomic sequence, given a chromosome, as start and an end position. So if the sequence of interest is the minus strand, one has to compute the reverse complement of the retrieved sequence, which can be done using functions provided in the matchprobes package.  The biomaRt vignette contains more examples on how to use this function.} 
\usage{getSequence( chromosome, start, end, id, type, seqType, upstream, downstream, mart, verbose=FALSE)}

\arguments{
  \item{chromosome}{Chromosome name}
  \item{start}{start position of sequence on chromosome}
  \item{end}{end position of sequence on chromosome}
  \item{id}{ An identifier or vector of identifiers.}
  \item{type}{The type of identifier used.  Supported types are hugo, ensembl, embl, entrezgene, refseq, ensemblTrans and unigene. Alternatively one can also use a filter to specify the type. Possible filters are given by the listFilters function}
  \item{seqType}{Type of sequence that you want to retrieve.  Allowed seqTypes are: cdna, peptide, 3utr, 5utr, genomic}
  \item{upstream}{To add the upstream sequence of a specified number of basepairs to the output.}
  \item{downstream}{To add the downstream sequence of a specified number of basepairs to the output.}
  \item{mart}{object of class Mart created using the useMart function}
  \item{verbose}{If verbose = TRUE then the XML query that was send to the webservice will be displayed.}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

seq = getSequence(id="BRCA1", type="hgnc_symbol", seqType="peptide", mart = mart)
show(seq)


seq = getSequence(id="1939_at", type="affy_hg_u95av2", seqType="gene_flank",upstream = 20, mart = mart)
show(seq)

}
}
\keyword{methods}

