\name{getSNP}
\alias{getSNP}
\title{Retrieves SNP's between a given chromosome start and end position}
\description{Retrieves SNP's given a chromosome name, a start and end position}
\usage{getSNP( chromosome, start, end, mart)}

\arguments{
  \item{chromosome}{Chromosome name}
  \item{start}{Start position on given chromosome}
  \item{end}{End position on given chromosome}
  \item{mart}{object of class Mart created using the useMart function}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
mart <- useMart("snp",dataset="hsapiens_snp")

snp = getSNP(chromosome = 8, start = 148350, end = 148612, mart = mart)
show(snp)

}
}
\keyword{methods}

