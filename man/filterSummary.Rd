\name{filterSummary}
\alias{filterSummary}
\title{Gives a summary of the filters groups and categories}
\description{Filters in BioMart databases are grouped together in filter groups and filter groups can in their turn form a higher collection in an filter category.  The filterSummary function gives a summary of the filter categories and groups present in the BioMart.  These group and category names can be used to display only a subset of the available filters in the listFilters function.}
\usage{filterSummary(mart)}

\arguments{
\item{mart}{object of class Mart, containing connections to the BioMart database.  You have to create such an object using the function useMart and update this Mart object by selecting a dataset.}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){

#marts <- listMarts()
#index<-grep("ensembl",marts)
#mart <- useMart(marts[index])
#datasets <- listDatasets(mart)
#mart <- useDataset(datasets[1],mart = mart)
#filterSummary(mart)
}
}
\keyword{methods}

