\name{useMart}
\alias{useMart}
\title{Connects to the selected BioMart}
\description{This function connects to the selected BioMart}
\usage{useMart(biomart, host, user, password, local = FALSE, mysql=FALSE)}

\arguments{
\item{biomart}{BioMart database name you want to connect to. Possible database names can be retrieved with the functio listMarts}
\item{mysql}{Boolean to specify if you want to access BioMart database via MySQL or via it's webservice of HTTP.  Default is using the webservice.}
\item{host}{If you want to use a local host or a miror database, use vector of hosts.  For connecting to public BioMarts this parameter is not required}
\item{user}{username, use vector of usernames for non public or miror BioMarts}
\item{password}{password, use vector of passwords for non-public or miror BioMarts}
\item{local}{boolean to specify if you want to use a locally installed BioMart or a mirror BioMart database}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){

#marts <- listMarts()
#index<-grep("ensembl",marts)
#mart <- useMart(marts[index])

#martDisconnect(mart = mart)
}
}
\keyword{methods}

