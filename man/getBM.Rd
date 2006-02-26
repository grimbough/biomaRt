\name{getBM}
\alias{getBM}
\title{Retrieves information from the BioMart database}
\description{This function retrieves information from the BioMart database}
\usage{getBM(attributes, filters, values, mart)}

\arguments{
  \item{attributes}{Attributes you want to retrieve.  A possible list of attributes can be retrieved using the function listAttributes.}
  \item{filters}{Filters that should be used in the query.  A possible list of filters can be retrieved using the function listFilters.}
  \item{values}{Values of the filter, e.g. list of affy id's}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function useMart and update the mart object by selecting a dataset using the function useDataset.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\examples{
if(interactive()){
mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart<-useDataset("hsapiens_gene_ensembl",mart)

getBM(attributes=c("external_gene_id","chr_name","band"),filters="affy_hg_u95av2",values=c("1939_at","1503_at","1454_at"), mart=mart)

martDisconnect(mart = mart)
}
}
\keyword{methods}

