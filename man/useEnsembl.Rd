\name{useEnsembl}
\alias{useEnsembl}
\alias{useEnsemblGenomes}
\title{Connects to the selected BioMart database and dataset hosted by Ensembl}
\description{A first step in using the biomaRt package is to select a BioMart
database and dataset to use.  The \code{useEnsembl} function enables one to 
connect to a specified BioMart database and dataset hosted by Ensembl without 
having to specify the Ensembl URL.  To know which BioMart databases are 
available see the \code{\link{listEnsembl}} and 
\code{\link{listEnsemblGenomes}} functions.  To know which datasets 
are available within a BioMart database, first select the BioMart database 
using \code{useEnsembl} and then use the \code{\link{listDatasets}} function 
on the selected Mart object.}

\usage{
useEnsembl(biomart, dataset, host, version = NULL, 
            GRCh = NULL, mirror = NULL, verbose = FALSE)
useEnsemblGenomes(biomart, dataset)
}

\arguments{
\item{biomart}{BioMart database name you want to connect to. Possible database names can be retrieved with the function \code{\link{listEnsembl}}}
\item{dataset}{Dataset you want to use.  To see the different datasets available within a biomaRt you can e.g. do:  mart = useEnsembl('genes'), followed by listDatasets(mart).}
\item{host}{Host to connect to.  Only needs to be specified if different from www.ensembl.org}
\item{version}{Ensembl version to connect to when wanting to connect to an archived Ensembl version}
\item{GRCh}{GRCh version to connect to if not the current GRCh38, currently this can only be 37}
\item{mirror}{Specify an Ensembl mirror to connect to.  The valid options here are 'www', 'uswest', 'useast', 'asia'.  If no mirror is specified the primary site at www.ensembl.org will be used.  Mirrors are not available for the Ensembl Genomes databases.}
\item{verbose}{Give detailed output of what the method is doing while in
  use, for debugging}
}

\author{Steffen Durinck & Mike Smith} 

\examples{
if(interactive()){

mart <- useEnsembl("ensembl")

## using the US West mirror
us_mart <- useEnsembl(biomart = "ensembl", mirror = "uswest")
}
}
\keyword{methods}

