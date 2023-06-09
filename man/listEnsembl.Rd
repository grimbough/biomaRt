\name{listEnsembl}
\alias{listEnsembl}
\alias{listEnsemblGenomes}
\title{lists the available BioMart databases hosted by Ensembl}
\description{This function returns a list of BioMart databases hosted by
Ensembl.  To establish a connection use the \code{\link{useEnsembl}} 
function.}
\usage{
    listEnsembl(mart = NULL, version = NULL, 
                GRCh = NULL, mirror = NULL, verbose = FALSE)
    listEnsemblGenomes(includeHosts = FALSE)
}
\arguments{
\item{mart}{mart object created with the useEnsembl function.  This is optional, as you usually use \code{\link{listMarts}} to see which marts there are to connect to.}
\item{version}{Ensembl version to connect to when wanting to connect to an archived Ensembl version}
\item{GRCh}{GRCh version to connect to if not the current GRCh38, currently this can only be 37}
\item{mirror}{Specify an Ensembl mirror to connect to.  The valid options here are 'www', 'useast', 'asia'.  If no mirror is specified the primary site at www.ensembl.org will be used.}
\item{verbose}{Give detailed output of what the method is doing, for debugging purposes}
\item{includeHosts}{If this option is set to \code{TRUE} a more detailed 
output is produced, including the URL used to access the corresponding mart.}
}

\author{Steffen Durinck, Mike L. Smith} 

\examples{
if(interactive()){
listEnsembl()
}
}
\keyword{methods}

