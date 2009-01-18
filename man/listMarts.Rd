\name{listMarts}
\alias{listMarts}
\title{lists the avilable BioMart databases}
\description{This function returns a list of BioMart databases to which biomaRt can connect to.  By default all public BioMart databases are displayed.  To establish a connection use the useMart function.}
\usage{listMarts(mart, host="www.biomart.org", path="/biomart/martservice", port=80, includeHosts = FALSE, archive=FALSE)}

\arguments{
\item{mart}{mart object created with the useMart function}
\item{host}{host to connect to if different then www.biomart.org}
\item{path}{path to martservice that should be pasted behind the host to get to web service URL}
\item{port}{port to use in HTTP communication} 
\item{includeHosts}{boolean to indicate if function should return host of the BioMart databases}
\item{archive}{Boolean to indicate if you want to access archived versions of BioMart database}
}
\author{Steffen Durinck, \url{http://www.stat.berkeley.edu/~steffen}} 

\seealso{}

\examples{
if(interactive()){
listMarts()
}
}
\keyword{methods}

