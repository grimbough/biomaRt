\name{listMarts}
\alias{listMarts}
\title{lists the avilable BioMart databases}
\description{This function returns a list of BioMart databases to which biomaRt can connect to.  By default all public BioMart databases are displayed.  To establish a connection use the useMart function.}
\usage{listMarts(mart, host="www.biomart.org", path="/biomart/martservice", port=80, includeHosts = FALSE, archive=FALSE, ssl.verifypeer=TRUE, verbose = FALSE)}

\arguments{
\item{mart}{mart object created with the useMart function}
\item{host}{host to connect to if different then www.biomart.org}
\item{path}{path to martservice that should be pasted behind the host to get to web service URL}
\item{port}{port to use in HTTP communication} 
\item{includeHosts}{boolean to indicate if function should return host of the BioMart databases}
\item{archive}{Boolean to indicate if you want to access archived versions of BioMart database}
\item{ssl.verifypeer}{Set SSL peer verification on or off.  By default ssl.verifypeer is set to TRUE}
\item{verbose}{Give detailed output of what the method is doing, for debugging purposes}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
listMarts()
}
}
\keyword{methods}

