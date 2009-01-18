\name{listFilters}
\alias{listFilters}
\title{lists the filters available in the selected dataset}
\description{Filters are what we use as inputs for a biomaRt query.  For example if we want to retrieve all entrezgene identifiers on chromosome X,  chromosome will be the filter with corresponding value X.  The listFilters functionlists the filters available in the selected dataset}
\usage{listFilters(mart, what = c("name","description"))}

\arguments{
\item{mart}{object of class Mart created using the useMart function}
\item{what}{vector which indicates what information to display about the available filters.  Acceptable values are name, description, fullDescription, type, page.}
}

\author{Steffen Durinck, http://www.stat.berkeley.edu/~steffen} 

\examples{

if(interactive()){
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
listFilters(mart)
}
}
\keyword{methods}

