\name{getBMlist}
\alias{getBMlist}
\title{Retrieves information from the BioMart database}
\description{This function is the main biomaRt query function.  Given a set of filters and corresponding values, it retrieves the user specified attributes from the BioMart database one is connected to}
\usage{getBMlist(attributes, filters = "", values = "", mart, list.names = NULL, na.value = NA, verbose = FALSE, giveWarning = TRUE)}

\arguments{
  \item{attributes}{Attributes you want to retrieve.  A possible list of attributes can be retrieved using the function listAttributes.}
  \item{filters}{Filters (one or more) that should be used in the query.  A possible list of filters can be retrieved using the function listFilters.}
  \item{values}{Values of the filter, e.g. vector of affy IDs.  If multiple filters are specified then the argument should be a list of vectors of which the position of each vector corresponds to the position of the filters in the filters argument.}
  \item{mart}{object of class Mart, created with the useMart function.}
\item{list.names}{names for objects in list}
\item{na.value}{value to give when result is empty}	 
\item{verbose}{When using biomaRt in webservice mode and setting verbose to TRUE, the XML query to the webservice will be printed.}
\item{giveWarning}{Gives a warning about best practices of biomaRt and recommends using getBM instead of getBMlist}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
mart <- useMart("ensembl")
datasets <- listDatasets(mart)

}
}
\keyword{methods}

