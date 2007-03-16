\name{attributeSummary}
\alias{attributeSummary}
\title{Gives a summary of the attributes groups and categories}
\description{Attributes in BioMart databases are grouped together in attribute groups and attribute groups can in their turn form a higher collection in an attribute category.  The attributeSummary function gives a summary of the attribute categories and groups present in the BioMart.  These group and category names can be used to display only a subset of the available attributes in the listAttributes function.}
\usage{attributeSummary(mart)}

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
#attributeSummary(mart)
}
}
\keyword{methods}

