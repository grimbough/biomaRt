\name{getPossibleXrefs}
\alias{getPossibleXrefs}
\title{Get possible cross-references in Ensembl}
\description{This function retrieves the possible cross-references present in Ensembl.  This is a very general function to see what can be extracted from Ensembl.  The results of this function can be used in the getXref function to extract the data of interest.}
\usage{getPossibleXrefs(mart)}

\arguments{
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Sean Davis, Steffen Durinck,  \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){
mart <- useMart("ensembl",mysql=TRUE)

xref <- getPossibleXrefs(mart = mart)
xref[1:10,]

martDisconnect(mart = mart)
}
}
\keyword{methods}

