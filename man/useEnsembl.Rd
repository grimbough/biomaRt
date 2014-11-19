\name{useEnsembl}
\alias{useEnsembl}
\title{Connects to the selected BioMart database and dataset hosted by Ensembl}
\description{A first step in using the biomaRt package is to select a BioMart database and dataset to use.  The useEnsembl function enables one to connect to a specified BioMart database and dataset hosted by Ensembl without having to specify the Ensembl URL.  To know which BioMart databases are available see the listEnsembl function.  To know which datasets are available within a BioMart database, first select the BioMart database using useMart and then use the listDatasets function on the selected BioMart, see listDatasets function.}
\usage{useEnsembl(biomart, dataset,host="www.ensembl.org", version = NULL, GRCh = NULL ,verbose = FALSE)}

\arguments{
\item{biomart}{BioMart database name you want to connect to. Possible database names can be retrieved with the functio listEnsembl}
\item{dataset}{Dataset you want to use.  To see the different datasets available within a biomaRt you can e.g. do:  mart = useEnsembl('ENSEMBL_MART_ENSEMBL'), followed by listDatasets(mart).}
\item{host}{Host to connect to if different then www.ensembl.org}
\item{version}{Ensembl version to connect to when wanting to connect to an archived Ensembl version}
\item{GRCh}{GRCh version to connect to if not the current GRCh38, currently this can only be 37}
\item{verbose}{Give detailed output of what the method is doing while in
  use, for debugging}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){

mart = useEnsembl("ENSEMBL_MART_ENSEMBL")
mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
}
}
\keyword{methods}

