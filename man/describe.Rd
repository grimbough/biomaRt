\name{describe}
\alias{describe}
\title{Retrieve a description of an attribute or filter}
\description{Retrieve a description of an attribute or filter}
\usage{describe(attribute,filter,mart)}

\arguments{
  \item{attribute}{An attribute name as shown by listAttributes}
  \item{filter}{A filter name as shown by listFilters}
  \item{mart}{object of class Mart, containing connections to the
    BioMart database.  You have to create such an object using the
    function useMart and update this Mart object by selecting a dataset.}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\examples{

if(interactive()){

mart <- useMart("ensembl","hsapiens_gene_ensembl")
describe(attribute="hgnc_symbol",mart=mart)
}
}
\keyword{methods}

