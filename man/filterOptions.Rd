\name{filterOptions}
\alias{filterOptions}
\title{Displays the filter options}
\description{Displays a set of predetermed values for the specified filter (if available).}
\usage{filterOptions(filter,mart)}

\arguments{
\item{filter}{A valid filter name.}
\item{mart}{object of class Mar created using the useMart functiont}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
filterOptions("chromosome_name", mart)
}
}
\keyword{methods}

