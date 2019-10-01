\name{listFilters}
\alias{listFilters}
\alias{searchFilters}
\title{List or search the filters available in the selected dataset}
\description{Filters are what we use as inputs for a biomaRt query.
  For example, if we want to retrieve all EntrezGene identifiers on
  chromosome X, \code{chromosome} will be the filter,
  with corresponding value X.}

\usage{
    listFilters(mart, what = c("name", "description"))
    searchFilters(mart, pattern)
}

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
\item{pattern}{Character vector defining the regular expression (\code{\link[base]{regex})}
to be used for the search.  If left blank the default is to use ".*" which will match everything.}
}

\author{Steffen Durinck, Mike Smith} 

\examples{

if(interactive()){

    ## list the available Ensembl marts and use Ensembl Genes
    listEnsembl()
    ensembl <- useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')

    ## list the available datasets in this Mart
    listFilters(mart = ensembl)

    ## the list of filters is long and not easy to read
    ## we can search for a term of interest to reduce this e.g. 'gene'
    searchFilters(mart = ensembl, pattern = "gene")
    
    ## search the available filters to find entries containing 'entrez' or 'hgnc'
    searchFilters(mart = ensembl, 'entrez|hgnc') 
}

}
\keyword{methods}

