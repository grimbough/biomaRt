\name{listFilters}
\alias{listFilters}
\title{lists the filters available in the selected dataset}
\description{lists the filters available in the selected dataset}
\usage{listFilters(mart, group, category, showGroups = FALSE)}

\arguments{
\item{mart}{object of class Mart, containing connections to the BioMart database.  You have to create such an object using the function useMart and update this Mart object by selecting a dataset.}
\item{group}{Show only the filters that belong to the specified filter group.  To get an overview of the filter groups set showGroups = TRUE or use the filterSummary function for a summary}
\item{category}{Show only the filters that belong to the specified filter category.  To get a summary of the filter category set showGroups = TRUE or use the filterSummary function for a summary}
\item{showGroups}{boolean to indicate if one wants to display the filter categories and groups along with their names and descriptions}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){
#marts <- listMarts()
#index<-grep("ensembl",marts)

#mart <- useMart(marts[index])
#datasets <- listDatasets(mart)
#mart <- useDataset(datasets[1],mart = mart)
#filters <- listFilters(mart)
#filters[1:10]
#martDisconnect(mart = mart)
}
}
\keyword{methods}

