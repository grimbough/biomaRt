\name{listDatasets}
\alias{listDatasets}
\alias{searchDatasets}

\title{List or search the datasets available in the selected BioMart database}
\description{Lists or search the datasets available in the selected BioMart database}

\usage{
    listDatasets(mart, verbose = FALSE)
    searchDatasets(mart, pattern)
}

\arguments{
    \item{mart}{object of class Mart created with the useMart function}
    \item{verbose}{Give detailed output of what the method is doing, for debugging purposes}
    \item{pattern}{Character vector defining the regular expression (\code{\link[base]{regex})}
    to be used for the search.  If left blank the default is to use ".*" which will match 
    everything and return the same as \code{listDatasets}.}
}

\author{Steffen Durinck, Mike Smith} 

\examples{

if(interactive()){

    ## list the available Ensembl marts and use Ensembl Genes
    listEnsembl()
    ensembl <- useEnsembl(biomart = "ensembl")

    ## list the available datasets in this Mart
    listDatasets(mart = ensembl)

    ## the list of Ensembl datasets grows ever larger (101 as of Ensembl 93)
    ## we can search for a term of interest to reduce the length e.g. 'sapiens'
    searchDatasets(mart = ensembl, pattern = "sapiens")
    
    ## search for any dataset containing the word Rat or rat
    searchDatasets(mart = ensembl, pattern = "(R|r)at")
}

}
\keyword{methods}

