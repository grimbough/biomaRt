\name{exportFASTA}
\alias{exportFASTA}
\title{Exports getSequence results to FASTA format}
\description{Exports getSequence results to FASTA format}
\usage{exportFASTA(martTable, file)}

\arguments{
\item{martTable}{A martTable that was the output of the getSequence function}
\item{file}{File to which you want to write the data}
 

}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{

if(interactive()){
mart <- useMart("ensembl")

#seq<-getSequence(chromosome=c(2,2),start=c(100000,30000),end=c(100300,30500),species="hsapiens",mart=mart)
#exportFASTA(seq,file="test.fasta")

martDisconnect(mart = mart)
}
}
\keyword{methods}

