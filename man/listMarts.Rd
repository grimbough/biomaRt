\name{listMarts}
\alias{listMarts}
\title{list the BioMart database covered by biomaRt package}
\description{This function returns a list for marts to which it can connect to. To establish a connection use the function martConnect with argument the name in the list (not the actual value e.g. use ensembl and not the full string of the database including release version).  You can also use this function to get the release information of the current databases which you using.}
\usage{listMarts(mart, host, user, password, port, includeHosts = FALSE, mysql=FALSE, archive=FALSE)}

\arguments{
\item{mart}{mart object}
\item{host}{host to connect to}
\item{user}{username}
\item{password}{password}
\item{port}{port to use in MySQL connection} 
\item{includeHosts}{boolean to indicate if function should return host of the BioMart databases}
\item{mysql}{Boolean to indicate to retrieve data from BioMart database via MySQL or via the BioMart webservice.  Default is via the webservice.}
\item{archive}{Boolean to indicate if you want to access archived versions od BioMart database
s.  This currently only works with mysql=TRUE and the ensembl BioMart database}
}
\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){
listMarts()
}
}
\keyword{methods}

