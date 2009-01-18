\name{useMart}
\alias{useMart}
\title{Connects to the selected BioMart database and dataset}
\description{A first step in using the biomaRt package is to select a BioMart database and dataset to use.  The useMart function enables one to connect to a specified BioMart database and dataset within this database.  To know which BioMart databases are available see the listMarts function.  To know which datasets are available within a BioMart database, first select the BioMart database using useMart and then use the listDatasets function on the selected BioMart, see listDatasets function.}
\usage{useMart(biomart, dataset, host="www.biomart.org", path="/biomart/martservice", port=80, archive=FALSE)}

\arguments{
\item{biomart}{BioMart database name you want to connect to. Possible database names can be retrieved with the functio listMarts}
\item{dataset}{Dataset you want to use.  To see the different datasets available within a biomaRt you can e.g. do:  mart = useMart('ensembl'), followed by listDatasets(mart).}
\item{host}{Host to connect to if different then www.biomart.org}
\item{path}{Path that should be pasted after to host to get access to the web service URL}
\item{port}{port to connect to, will be pasted between host and path}
\item{archive}{Boolean to indicate if you want to access archived versions od BioMart databases.  This currently only works with mysql=TRUE and the ensembl BioMart database}
}

\author{Steffen Durinck,\url{http://www.stat.berkeley.edu/~steffen}} 

\examples{
if(interactive()){

mart = useMart("ensembl")
mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
}
}
\keyword{methods}

