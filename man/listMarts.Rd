\name{listMarts}
\alias{listMarts}
\title{lists the available BioMart databases}
\description{This function returns a list of BioMart databases to which 
biomaRt can connect.  By default the Ensembl BioMart databases are displayed.
To establish a connection use the \link{useMart} function.}

\usage{listMarts(mart = NULL, host="https://www.ensembl.org", path="/biomart/martservice",
port, includeHosts = FALSE, archive = FALSE, httr_config, verbose = FALSE)}

\arguments{
\item{mart}{mart object created with the \code{\link{useMart}} function.  
This is optional, as you usually use \code{\link{listMarts}} to see which 
marts there are to connect to.}
\item{host}{Host to connect to. Defaults to \code{www.ensembl.org}}
\item{path}{path to martservice that should be pasted behind the host to get to web service URL}
\item{port}{port to use in HTTP communication} 
\item{includeHosts}{boolean to indicate if function should return host of the BioMart databases}
\item{archive}{Boolean to indicate if you want to access archived versions of 
BioMart database. Note that this argument is now defunct and setting this value
to \code{TRUE} will produce an error.
A better alternative is to specify the url of the archived 
BioMart you want to access.  For Ensembl you can view the list of archives 
using \code{\link{listEnsemblArchives}}}
\item{httr_config}{Some hosts require specific HTTP settings to be used when 
connecting. This argument takes the output of \code{\link[httr]{config}} and
will be used when connecting to \code{host}.  Can be ignored if you experience
no problems accessing \code{host}.}
\item{verbose}{Give detailed output of what the method is doing, for debugging purposes.}
}

\details{
If you receive an error message saying 'Unexpected format to the list of 
available marts', this is often because there is a problem with the BioMart
server you are trying to connect to, and something other than the list of 
available marts is being returned - often some like a 'down for maintainance'
page.  If you browse to the provided URL and find a page that starts with
'\code{<MartRegistry>}' this is the correct listing and you should report
the issue on the Bioconductor support site: https://support.bioconductor.org
}

\author{Steffen Durinck, Mike Smith} 

\examples{
if(interactive()){
listMarts()
}
}
\keyword{methods}

