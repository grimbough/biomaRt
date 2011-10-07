\name{filterType}
\alias{filterType}
\title{Displays the filter type}
\description{Displays the type of the filer given a filter name.}
\usage{filterType(filter,mart)}

\arguments{
\item{filter}{A valid filter name. Valid filters are given by the listFilters function}
\item{mart}{object of class Mart, created using the useMart function}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
filterType("chromosome_name", mart)
}
}
\keyword{methods}

