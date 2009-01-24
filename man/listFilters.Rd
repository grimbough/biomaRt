\name{listFilters}
\alias{listFilters}
\title{lists the filters available in the selected dataset}
\description{Filters are what we use as inputs for a biomaRt query.
  For example, if we want to retrieve all entrezgene identifiers on
  chromosome X, \code{chromosome} will be the filter,
  with corresponding value X.}

\usage{listFilters(mart, what = c("name", "description"),
                   group = "DEFUNCT")}

\arguments{
\item{mart}{object of class \code{Mart} created using the \code{\link{useMart}} function}
\item{what}{character vector indicating what information to display about the
  available filters.  Valid values are
  \code{name},
  \code{description},
  \code{options},
  \code{fullDescription},
  \code{filters},
  \code{type},
  \code{operation},
  \code{filters8},
  \code{filters9}.}
\item{group}{defunct. If you need advice how to replace that
  functionality, please contact the package maintainer for advice.}
}

\author{Steffen Durinck, http://www.stat.berkeley.edu/~steffen} 

\examples{

if(interactive()){
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
listFilters(mart)
}

}
\keyword{methods}

