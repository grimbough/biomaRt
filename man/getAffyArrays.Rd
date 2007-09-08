\name{getAffyArrays}
\alias{getAffyArrays}
\title{Displays affy identifiers of Affymetrix chips supported by the Ensembl package}
\description{Displays affy identifiers of Affymetrix chips supported by the Ensembl package}
\usage{getAffyArrays( mart )}

\arguments{
  \item{mart}{object of class Mart, representing connections to the
    BioMart databases}
}

\author{Steffen Durinck, \url{http://www.stat.berkeley.edu/~steffen}} 

\examples{
if(interactive()){
mart <- useMart("ensembl",mysql=TRUE)
getAffyArrays(mart = mart)
martDisconnect( mart )
}
}
\keyword{methods}

