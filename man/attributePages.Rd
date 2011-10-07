\name{attributePages}
\alias{attributePages}
\title{Gives a summary of the attribute pages}
\description{Attributes in BioMart databases are grouped together in attribute pages.  The attributePages function gives a summary of the attribute categories and groups present in the BioMart.  These page names can be used to display only a subset of the available attributes in the listAttributes function.}
\usage{attributePages(mart)}

\arguments{
\item{mart}{object of class Mart, created with the useMart function.}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
attributeSummary(mart)
}
}
\keyword{methods}

