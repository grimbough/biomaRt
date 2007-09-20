\name{listMarts}
\alias{listMarts}
\title{lists the avilable BioMart databases}
\description{This function returns a list of BioMart databases to which biomaRt can connect to.  By default all public BioMart databases are displayed.  To establish a connection use the useMart function.}
\usage{listMarts(mart, host, user, password, port, includeHosts = FALSE, mysql=FALSE, archive=FALSE)}

\arguments{
\item{mart}{mart object created with the useMart function}
\item{host}{host to connect to}
\item{user}{username}
\item{password}{password}
\item{port}{port to use in MySQL connection} 
\item{includeHosts}{boolean to indicate if function should return host of the BioMart databases}
\item{mysql}{Boolean to indicate to retrieve data from BioMart database via MySQL or via the BioMart webservice.  Default is via the webservice.}
\item{archive}{Boolean to indicate if you want to access archived versions od BioMart database
s.  This currently only works with mysql=TRUE and the ensembl BioMart database}
}
\author{Steffen Durinck, \url{http://www.stat.berkeley.edu/~steffen}} 

\seealso{}

\examples{
if(interactive()){
listMarts()
}
}
\keyword{methods}

