\name{getLDS}
\alias{getLDS}
\title{Retrieves information from two linked datasets}
\description{This function retrieves information from two linked BioMart datasets}
\usage{getLDS(attributes, filters = "", values = "", mart, attributesL, filtersL = "", valuesL = "", martL, verbose = FALSE, uniqueRows = TRUE)}
\arguments{
  \item{attributes}{Attributes you want to retrieve of primary dataset.  A possible list of attributes can be retrieved using the function listAttributes.}
  \item{filters}{Filters that should be used in the query. These filters will be applied to primary dataset.  A possible list of filters can be retrieved using the function listFilters.}
  \item{values}{Values of the filter, e.g. list of affy IDs}
  \item{mart}{object of class Mart, representing the primary dataset on which filters apply .You can create such an object using the function useMart and update the mart object by selecting a dataset using the function useDataset.}
  \item{attributesL}{Attributes of linked dataset that needs to be retrieved}
  \item{filtersL}{Filters to be applied to the linked dataset}
  \item{valuesL}{Values for the linked dataset filters}
  \item{martL}{Mart object representing linked dataset}
 \item{verbose}{When using biomaRt in webservice mode and setting verbose to TRUE, the XML query to the webservice will be printed.  Alternatively in MySQL mode the MySQL query will be printed.}
  \item{uniqueRows}{Logical to indicate if the BioMart web service should return unique rows only or not.  Has the value of either TRUE or FALSE}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl") 
getLDS(attributes = c("hgnc_symbol","chromosome_name", "start_position"), filters = "hgnc_symbol", values = "TP53", mart = human, attributesL = c("chromosome_name","start_position"), martL = mouse)
}
}
\keyword{methods}

