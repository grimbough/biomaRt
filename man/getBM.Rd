\name{getBM}
\alias{getBM}
\title{Retrieves information from the BioMart database}
\description{This function retrieves information from the BioMart database}
\usage{getBM(attributes, filters, values, mart, curl = NULL, output = "data.frame", list.names = NULL, na.value = NA)}

\arguments{
  \item{attributes}{Attributes you want to retrieve.  A possible list of attributes can be retrieved using the function listAttributes.}
  \item{filters}{Filters that should be used in the query.  A possible list of filters can be retrieved using the function listFilters.}
  \item{values}{Values of the filter, e.g. list of affy IDs}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function useMart and update the mart object by selecting a dataset using the function useDataset.}
  \item{curl}{An optional 'CURLHandle' object, that can be used to speed up getBM when used in a loop.}
  \item{output}{Determines the output of getBM which can be either a data.frame (default) or a list.}
  \item{list.names}{In case a list was selected as output the different elements in the list can be given user defined names with this argument}
  \item{na.value}{In case a list was selected as output, the value of na.value will used when missing elements are present in the output}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\examples{
if(interactive()){
mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart<-useDataset("hsapiens_gene_ensembl",mart)
getBM(attributes=c("affy_hg_u95av2","hgnc_symbol","chromosome_name","band"),filters="affy_hg_u95av2",values=c("1939_at","1503_at","1454_at"), mart=mart)

}
}
\keyword{methods}

