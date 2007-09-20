\name{listFilters}
\alias{listFilters}
\title{lists the filters available in the selected dataset}
\description{Filters are what we use as inputs for a biomaRt query.  For example if we want to retrieve all entrezgene identifiers on chromosome X,  chromosome will be the filter with corresponding value X.  The listFilters functionlists the filters available in the selected dataset}
\usage{listFilters(mart, group, category, showGroups = FALSE, showType = FALSE)}

\arguments{
\item{mart}{object of class Mart created using the useMart function}
\item{group}{Show only the filters that belong to the specified filter group.  To get an overview of the filter groups set showGroups = TRUE or use the filterSummary function for a summary}
\item{category}{Show only the filters that belong to the specified filter category.  To get a summary of the filter category set showGroups = TRUE or use the filterSummary function for a summary}
\item{showGroups}{Boolean to indicate if one wants to display the filter categories and groups along with their names and descriptions}
\item{showType}{Boolean to indicate if one wants to display the type of filter (boolean, text,list) along with their names and descriptions.  This is especially useful to see wether a filter is a boolean or not.}
}

\author{Steffen Durinck, http://www.stat.berkeley.edu/~steffen} 

\examples{

if(interactive()){
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
listFilters(mart)
}
}
\keyword{methods}

