\name{listMarts}
\alias{listMarts}
\title{list the BioMart database covered by biomaRt package}
\description{This function returns a list for marts to which it can connect to. To establish a connection use the function martConnect with argument the name in the list (not the actual value e.g. use ensembl and not the full string of the database including release version).  You can also use this function to get the release information of the current databases which you using.}
\usage{listMarts()}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
listMarts()
}
\keyword{methods}

