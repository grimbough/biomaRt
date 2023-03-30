\name{getBM}
\alias{getBM}
\title{Retrieves information from the BioMart database}
\description{This function is the main biomaRt query function.  Given a set of filters and corresponding values, it retrieves the user specified attributes from the BioMart database one is connected to.}
\usage{getBM(attributes, filters = "", values = "", mart, curl = NULL, 
checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE,
quote = "\"", useCache = TRUE)}
\value{A \code{data.frame}. There is no implicit mapping between its rows and the function arguments (e.g. \code{filters}, \code{values}), therefore make sure to have the relevant identifier(s) returned by specifying them in \code{attributes}. See Examples.}

\arguments{
  \item{attributes}{Attributes you want to retrieve.  A possible list of attributes can be retrieved using the function listAttributes.}
  \item{filters}{Filters (one or more) that should be used in the query.  A possible list of filters can be retrieved using the function listFilters.}
  \item{values}{Values of the filter, e.g. vector of affy IDs.  If multiple filters are specified then the argument should be a list of vectors of which the position of each vector corresponds to the position of the filters in the filters argument.}
  \item{mart}{object of class Mart, created with the useMart function.}
  \item{curl}{An optional 'CURLHandle' object, that can be used to speed up getBM when used in a loop.}
  \item{checkFilters}{Sometimes attributes where a value needs to be specified, for example upstream_flank with value 20 for obtaining upstream sequence flank regions of length 20bp, are treated as filters in BioMarts.  To enable such a query to work, one must specify the attribute as a filter and set checkFilters = FALSE for the query to work.}
 \item{verbose}{When using biomaRt in webservice mode and setting verbose to TRUE, the XML query to the webservice will be printed.}
 \item{uniqueRows}{If the result of a query contains multiple identical rows, setting this argument to TRUE (default) will result in deleting the duplicated rows in the query result at the server side.}
 \item{bmHeader}{Boolean to indicate if the result retrieved from the
  BioMart server should include the data headers or not, defaults to FALSE.  This should
  only be switched on if the default behavior results in errors,
  setting to on might still be able to retrieve your data in that case}
 \item{quote}{Sometimes parsing of the results fails due to errors in the Ensembl data fields such as containing a quote, in such cases you can try to change the value of quote to try to still parse the results.}
 \item{useCache}{Boolean indicating whether the results cache should be used.  
 Setting to \code{FALSE} will disable reading and writing of the cache.  This
 argument is likely to disappear after the cache functionality has been tested
 more thoroughly.}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
  mart <- useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl")
                     
  getBM(attributes = c("affy_hg_u95av2", "hgnc_symbol", "chromosome_name", "band"),
        filters    = "affy_hg_u95av2",
        values     = c("1939_at","1503_at","1454_at"), 
        mart       = mart)
  }
}
\keyword{methods}

