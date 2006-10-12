\name{exportFASTA}
\alias{exportFASTA}
\title{Exports getSequence results to FASTA format}
\description{Exports getSequence results to FASTA format}
\usage{exportFASTA(sequences, file)}

\arguments{
\item{sequences}{A data.frame that was the output of the getSequence function}
\item{file}{File to which you want to write the data}
 

}

\author{Steffen Durinck} 

\examples{

if(interactive()){
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#seq<-getSequence(chromosome=c(2,2),start=c(100000,30000),end=c(100300,30500),mart=mart)
#exportFASTA(seq,file="test.fasta")

martDisconnect(mart = mart)
}
}
\keyword{methods}

