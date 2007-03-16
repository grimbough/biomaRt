\name{listAttributes}
\alias{listAttributes}
\title{lists the attributes available in the selected dataset}
\description{lists the datasets available in the selected dataset}
\usage{listAttributes(mart, group, category, showGroups = FALSE)}

\arguments{
\item{mart}{object of class Mart, containing connections to the BioMart database.  You have to create such an object using the function useMart and update this Mart object by selecting a dataset.}
\item{group}{Show only the attributes that belong to the specified attribute group.  To get a summary of the attribute groups set showGroups = TRUE or use the attributeSummary function}
\item{category}{Show only the attributes that belong to the specified attribute category.  To get a summary of the attribute pages set showGroups = TRUE or use the attributeSummary function}
\item{showGroups}{boolean to indicate if one wants to display the attribute categories and groups along with their names and descriptions}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){

#marts <- listMarts()
#index<-grep("ensembl",marts)
#mart <- useMart(marts[index])
#datasets <- listDatasets(mart)
#mart <- useDataset(datasets[1],mart = mart)
#attributes <- listAttributes(mart)
#attributes[1:10]
#martDisconnect(mart = mart)
}
}
\keyword{methods}

