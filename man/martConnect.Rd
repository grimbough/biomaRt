\name{martConnect}
\alias{martConnect}
\title{Connects to available BioMarts}
\description{This function connects to available BioMarts}
\usage{martConnect(biomarts = "ensembl", host = NULL, user = NULL, password = NULL, mart = NULL, local = FALSE)}

\arguments{
\item{biomarts}{vector with BioMart database names you want to connect to. Use one or a vector of the following BioMart database names:  ensembl (default, for gene information, GO and OMIM annotation),sequence (for sequence retrieval), snp (for SNP retrieval), uniprot (for protein related data), vega (for VEGA gene annotation)}
\item{host}{If you want to use a local host or a miror database, use vector of hosts.  For connecting to public BioMarts this parameter is not required}
\item{user}{username, use vector of usernames for non public or miror BioMarts}
\item{password}{password, use vector of passwords for non-public or miror BioMarts}
\item{mart}{Mart object, use if you want to add a connection to a BioMart to an existing Mart object}
\item{local}{boolean to specify if you want to use a locally installed BioMart or a mirror BioMart database}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){
mart <- martConnect(biomarts = c("ensembl","sequence"))

martDisconnect(mart = mart)
}
}
\keyword{methods}

