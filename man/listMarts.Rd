\name{listMarts}
\alias{listMarts}
\title{lists the available BioMart databases}
\description{This function returns a list of BioMart databases to which biomaRt can connect to.  By default all public BioMart databases are displayed.  To establish a connection use the useMart function.}
\usage{listMarts(mart = NULL, host="www.ensembl.org", path="/biomart/martservice",
port=80, includeHosts = FALSE, archive = FALSE, ssl.verifypeer = TRUE, 
ensemblRedirect = TRUE, verbose = FALSE)}

\arguments{
\item{mart}{mart object created with the \code{\link{useMart}} function.  
This is optional, as you usually use \code{\link{listMarts}} to see which 
marts there are to connect to.}
\item{host}{Host to connect to. Defaults to \code{www.ensembl.org}}
\item{path}{path to martservice that should be pasted behind the host to get to web service URL}
\item{port}{port to use in HTTP communication} 
\item{includeHosts}{boolean to indicate if function should return host of the BioMart databases}
\item{archive}{Boolean to indicate if you want to access archived versions of BioMart database}
\item{ssl.verifypeer}{Set SSL peer verification on or off.  By default ssl.verifypeer is set to TRUE}
\item{ensemblRedirect}{By default when you access Ensembl BioMart it will 
redirect you to your local mirror, even if you have set a region specific 
mirror in the \code{host} argument.  By setting this argument to \code{FALSE} 
you can override this behaviour and force access to your specified \code{host}.
Defaults to \code{TRUE} to prevent overwhelming the main Ensembl site.
If you are accessing a BioMart instance other than Ensembl this should have
no effect.}
\item{verbose}{Give detailed output of what the method is doing, for debugging purposes}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
listMarts()
}
}
\keyword{methods}

