\name{get}
\alias{get}
\title{Retrieves information from the BioMart database}
\description{This function retrieves information from the BioMart database}
\usage{get(attributes, filter, values, mart)}

\arguments{
  \item{attributes}{Attributes you want to retrieve.  A possible list of attributes can be retrieved using the function listAttributes.}
  \item{filter}{Filter that should be used in the query.  A possible list of filters can be retrieved using the function listFilters.}
  \item{values}{Values of the filter, e.g. list of affy id's}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function useMart and update the mart object by selecting a dataset using the function useDataset.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{
if(interactive()){
#marts<-listMarts()
#index<-grep("ensembl",marts)
#mart <- useMart(marts[index])
#datasets <- listDatasets(mart)
#mart<-useDataset("hsapiens_gene_ensembl",mart)

#get(attributes=c("external_gene_id","chr_name","band"),filter="affy_hg_u95_v2",values=c("1939_at","1503_at","1454_at"), mart=mart)

#martDisconnect(mart = mart)
}
}
\keyword{methods}

