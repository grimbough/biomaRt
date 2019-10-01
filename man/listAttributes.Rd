\name{listAttributes}
\alias{listAttributes}
\alias{searchAttributes}
\title{lists the attributes available in the selected dataset}
\description{Attributes are the outputs of a biomaRt query, they are the 
information we want to retrieve.  For example if we want to retrieve all 
EntrezGene identifiers of genes located on chromosome X, \code{entrezgene_id}
will be the attribute we use in the query.  The \code{listAttributes} function 
lists the available attributes in the selected dataset.}

\usage{
    listAttributes(mart, page,what = c("name","description","page"))
    searchAttributes(mart, pattern)
}

\arguments{
    \item{mart}{object of class Mart created using the useMart function}
    \item{page}{Show only the attributes that belong to the specified attribute page.}
    \item{what}{vector of types of information about the attributes that need to be displayed.  
    Can have values like name, description, fullDescription, page}
    \item{pattern}{Character vector defining the regular expression (\code{\link[base]{regex})}
    to be used for the search.  If left blank the default is to use ".*" which will match 
    everything.}
}

\author{Steffen Durinck, Mike Smith} 

\examples{

if(interactive()){

    ## list the available Ensembl marts and use Ensembl Genes
    listEnsembl()
    ensembl <- useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')

    ## list the available datasets in this Mart
    listAttributes(mart = ensembl)

    ## the list of attributes is very long and gets truncated by R
    ## we can search for a term of interest to filter this e.g. 'start'
    searchAttributes(mart = ensembl, pattern = "start")
    
    ## filter the attributes to give only entries containing 'entrez' or 'hgnc'
    searchAttributes(mart = ensembl, 'entrez|hgnc') 
}

}
\keyword{methods}

