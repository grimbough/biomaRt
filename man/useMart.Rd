\name{useMart}
\alias{useMart}
\title{Connects to the selected BioMart database and dataset}
\description{A first step in using the biomaRt package is to select a BioMart database and dataset to use.  The useMart function enables one to connect to a specified BioMart database and dataset within this database.  To know which BioMart databases are available see the listMarts function.  To know which datasets are available within a BioMart database, first select the BioMart database using useMart and then use the listDatasets function on the selected BioMart, see listDatasets function.}
\usage{useMart(biomart, dataset, host="https://www.ensembl.org",
path="/biomart/martservice", port, archive=FALSE, version, verbose = FALSE)}

\arguments{
\item{biomart}{BioMart database name you want to connect to. Possible database names can be retrieved with the functio listMarts}
\item{dataset}{Dataset you want to use.  To see the different datasets available within a biomaRt you can e.g. do:  mart = useMart('ensembl'), followed by listDatasets(mart).}
\item{host}{Host to connect to. Defaults to \code{www.ensembl.org}}
\item{path}{Path that should be pasted after to host to get access to the web service URL}
\item{port}{port to connect to, will be pasted between host and path}
\item{archive}{Boolean to indicate if you want to access archived versions of BioMart databases.  Note that this argument is now deprecated and will be removed in the future.  A better alternative is to leave archive = FALSE and to specify the url of the archived BioMart you want to access.  For Ensembl you can view the list of archives using \code{\link{listEnsemblArchives}}}
\item{version}{Use version name instead of biomart name to specify which BioMart you want to use}
\item{verbose}{Give detailed output of what the method is doing while in
  use, for debugging}
}

\author{Steffen Durinck, Mike L. Smith} 

\examples{
if(interactive()){

    mart = useMart("ensembl")
    mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
}
}
\keyword{methods}

