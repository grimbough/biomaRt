\name{martDisconnect}
\alias{martDisconnect}
\title{Close a connection with to BioMart database}
\description{This function closes a connection with a BioMart database and returns TRUE/FALSE upon completion}
\usage{martDisconnect(mart = NULL)}

\arguments{
  \item{mart}{object of class Mart, representing a connection to a BioMart database}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
#connect to Ensembl
mart <- martConnect("ensembl")

#and disconnect ....
martDisconnect(mart)

}
\keyword{methods}

