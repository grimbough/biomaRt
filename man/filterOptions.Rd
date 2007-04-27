\name{filterOptions}
\alias{filterOptions}
\title{Displays the filter options}
\description{Displays a set of predetermed values for the specified filter (if available).}
\usage{filterOptions(filter,mart)}

\arguments{
\item{filter}{A valid filter name.}
\item{mart}{object of class Mart, containing connections to the BioMart database.  You have to create such an object using the function useMart and update this Mart object by selecting a dataset.}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){
#mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
#filterOptions("chromosome_name", mart)
}
}
\keyword{methods}

