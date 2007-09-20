\name{filterSummary}
\alias{filterSummary}
\title{Gives a summary of the filters groups and categories}
\description{Filters in BioMart databases are grouped together in filter groups and filter groups can in their turn form a higher collection in an filter category.  The filterSummary function gives a summary of the filter categories and groups present in the BioMart.  These group and category names can be used to display only a subset of the available filters in the listFilters function.}
\usage{filterSummary(mart)}

\arguments{
\item{mart}{object of class Mart created with the useMart function.}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){
mart=useMart("ensembl", dataset="hsapiens_gene_ensembl")
filterSummary(mart)
}
}
\keyword{methods}

