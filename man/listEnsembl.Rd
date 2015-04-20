\name{listEnsembl}
\alias{listEnsembl}
\title{lists the available BioMart databases hosted by Ensembl}
\description{This function returns a list of BioMart databases hosted by Ensembl.  To establish a connection use the useMart function.}
\usage{listEnsembl(mart = NULL,host = "www.ensembl.org", version = NULL, GRCh = NULL, mirror = NULL, verbose = FALSE)}

\arguments{
\item{mart}{mart object created with the useEnsembl function.  This is optional, as you usually use listMarts to see which marts there are to connect to.}
\item{host}{Host to connect to if different then www.ensembl.org}
\item{version}{Ensembl version to connect to when wanting to connect to an archived Ensembl version}
\item{GRCh}{GRCh version to connect to if not the current GRCh38, currently this can only be 37}
\item{mirror}{Mirror to connect to [uswest,useast,asia], If no mirror is specified the default will be used which is www.ensembl.org hosted in the UK}
\item{verbose}{Give detailed output of what the method is doing, for debugging purposes}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
listEnsembl()
}
}
\keyword{methods}

