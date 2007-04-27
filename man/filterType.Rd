\name{filterType}
\alias{filterType}
\title{Displays the filter type}
\description{Displays the type of the filer given a filter name.  Is useful to see which filters are boolean and where TRUE and FALSE is needed as values.}
\usage{filterType(filter,mart)}

\arguments{
\item{filter}{A valid filter name.}
\item{mart}{object of class Mart, containing connections to the BioMart database.  You have to create such an object using the function useMart and update this Mart object by selecting a dataset.}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){
#mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
#filterType("chromosome_name", mart)
}
}
\keyword{methods}

