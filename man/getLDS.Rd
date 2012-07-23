\name{getLDS}
\alias{getLDS}
\title{Retrieves information from two linked datasets}
\description{This function is the main biomaRt query function that links 2 datasets and retrieves information from these linked BioMart datasets.  In Ensembl this translates to homology mapping.}
\usage{getLDS(attributes, filters = "", values = "", mart, attributesL,
filtersL = "", valuesL = "", martL, verbose = FALSE, uniqueRows = TRUE, bmHeader=TRUE)}
\arguments{
  \item{attributes}{Attributes you want to retrieve of primary dataset.  A possible list of attributes can be retrieved using the function listAttributes.}
  \item{filters}{Filters that should be used in the query. These filters will be applied to primary dataset.  A possible list of filters can be retrieved using the function listFilters.}
  \item{values}{Values of the filter, e.g. list of affy IDs}
  \item{mart}{object of class Mart created with the useMart function.}
  \item{attributesL}{Attributes of linked dataset that needs to be retrieved}
  \item{filtersL}{Filters to be applied to the linked dataset}
  \item{valuesL}{Values for the linked dataset filters}
  \item{martL}{Mart object representing linked dataset}
 \item{verbose}{When using biomaRt in webservice mode and setting verbose to TRUE, the XML query to the webservice will be printed.  Alternatively in MySQL mode the MySQL query will be printed.}
  \item{uniqueRows}{Logical to indicate if the BioMart web service should return unique rows only or not.  Has the value of either TRUE or FALSE}
  \item{bmHeader}{Boolean to indicate if the result retrieved from the
      BioMart server should include the data headers or not, defaults to
      TRUE.  This should only be switched off if the default behavior
      results in errors, setting to off might still be able to retrieve
      your data in that case}
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

