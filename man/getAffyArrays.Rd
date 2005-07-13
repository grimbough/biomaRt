\name{getAffyArrays}
\alias{getAffyArrays}
\title{Displays affy identifiers of Affymetrix chips supported by the Ensembl package}
\description{Displays affy identifiers of Affymetrix chips supported by the Ensembl package}
\usage{getAffyArrays( mart )}

\arguments{
  \item{mart}{object of class Mart, representing connections to the
    BioMart databases}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){
mart <- martConnect()
getAffyArrays(mart = mart)
martDisconnect( mart )
}
}
\keyword{methods}

