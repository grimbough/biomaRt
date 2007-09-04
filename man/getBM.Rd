\name{getBM}
\alias{getBM}
\title{Retrieves information from the BioMart database}
\description{This function retrieves information from the BioMart database}
\usage{getBM(attributes, filters = "", values = "", mart, curl = NULL, output = "data.frame", list.names = NULL, na.value = NA, checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE)}

\arguments{
  \item{attributes}{Attributes you want to retrieve.  A possible list of attributes can be retrieved using the function listAttributes.}
  \item{filters}{Filters (one or more) that should be used in the query.  A possible list of filters can be retrieved using the function listFilters.}
  \item{values}{Values of the filter, e.g. vector of affy IDs.  If multiple filters are specified then the argument should be a list of vectors of which the position of each vector corresponds to the position of the filters in the filters argument.}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can creat such an object using the function useMart and update the mart object by selecting a dataset using the function useDataset.}
  \item{curl}{An optional 'CURLHandle' object, that can be used to speed up getBM when used in a loop.}
  \item{output}{Determines the output of getBM which can be either a data.frame (default) or a list.}
  \item{list.names}{In case a list was selected as output the different elements in the list can be given user defined names with this argument}
  \item{na.value}{In case a list was selected as output, the value of na.value will used when missing elements are present in the output}
  \item{checkFilters}{Sometimes attributes where a value needs to be specified, for example upstream\_flank with value 20 for obtaining upstream sequence flank regions of length 20bp, are treated as filters in BioMarts.  To enable such a query to work, one must specify the attribute as a filter and set checkFilters = FALSE for the query to work.}
 \item{verbose}{When using biomaRt in webservice mode and setting verbose to TRUE, the XML query to the webservice will be printed.  Alternatively in MySQL mode the MySQL query will be printed.}
 \item{uniqueRows}{If the result of a query contains multiple identical rows, setting this argument to TRUE (default) will result in deleting the duplicated rows in the query result at the server side.}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart<-useDataset("hsapiens_gene_ensembl",mart)
getBM(attributes=c("affy_hg_u95av2","hgnc_symbol","chromosome_name","band"),filters="affy_hg_u95av2",values=c("1939_at","1503_at","1454_at"), mart=mart)

}
}
\keyword{methods}

