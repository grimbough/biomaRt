\name{attributeSummary}
\alias{attributeSummary}
\title{Gives a summary of the attributes groups and categories}
\description{Attributes in BioMart databases are grouped together in attribute groups and attribute groups can in their turn form a higher collection in an attribute category.  The attributeSummary function gives a summary of the attribute categories and groups present in the BioMart.  These group and category names can be used to display only a subset of the available attributes in the listAttributes function.}
\usage{attributeSummary(mart)}

\arguments{
\item{mart}{object of class Mart, created with the useMart function.}
}

\author{Steffen Durinck, http://www.stat.berkeley.edu/~steffen} 

\examples{

if(interactive()){
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
attributeSummary(mart)
}
}
\keyword{methods}

