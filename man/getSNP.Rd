\name{getSNP}
\alias{getSNP}
\title{Retrieves SNP's between a given chromosome start and end position}
\description{Retrieves SNP's between a given chromosome start and end position}
\usage{getSNP( chromosome, start, end, mart)}

\arguments{
  \item{chromosome}{Chromosome name}
  \item{start}{start position on given chromosome}
  \item{end}{end position on given chromosome}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
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

