\name{listFilterOptions}
\alias{searchFilterOptions}
\alias{listFilterOptions}
\title{List or search the options available for a specified filter.}
\description{Some filters have a predefined list of values that can be used
to search them.  These functions give access to this list of options for
a named filter, so you can check in the case where your biomaRt query is 
not finding anything.}

\usage{
    listFilterOptions(mart, filter)
    searchFilterOptions(mart, filter, pattern)
}

\arguments{
    \item{mart}{object of class \code{Mart} created using the \code{\link{useMart}}, 
    or \code{\link{useEnsembl}} functions}
    \item{filter}{The name of the filter whose options should be listed or searched.  You
    can list available filters via \code{\link{listFilters}}}
    \item{pattern}{Character vector defining the regular expression (\code{\link[base]{regex})}
    to be used for the search.  If left blank the default is to use ".*" which will match everything.}
}

\author{Mike Smith} 

\examples{

if(interactive()){

    ## Use the Ensembl human genes dataset
    ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

    ## we can search for the name of a filter we're interested in e.g. 'phenotype'
    ## we need to use the name of the filter in the next function
    searchFilters(ensembl, pattern = "phenotype")
    
    ## list all the options available to the 'phenotype_source' filter
    listFilterOptions(mart = ensembl, filter = "phenotype_source")
    
    ## search the 'phenotype_description' filter for the term 'crohn'
    searchFilterOptions(mart = ensembl, 
                       filter = "phenotype_description", 
                       pattern = "crohn")
}

}

\seealso{
    \code{\link{listFilters}}
}

\keyword{methods}

