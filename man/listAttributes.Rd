\name{listAttributes}
\alias{listAttributes}
\title{lists the attributes available in the selected dataset}
\description{Attributes are the outputs of a biomaRt query, they are the information we want to retrieve.  For example if we want to retrieve all entrez gene identifiers of genes located on chromosome X, entrezgene will be the attribute we use in the query.  The listAttributes function lists the available attributes in the selected dataset}
\usage{listAttributes(mart, page,what = c("name","description"), group, category, showGroups = FALSE)}

\arguments{
\item{mart}{object of class Mart created using the useMart function}
\item{page}{Show only the attributes that belong to the specified attribute page.}
\item{what}{vector of types of information about the attributes that need to be displayed.  Can have values like name, description, fullDescription, page}
\item{group}{Availability of group argument is pending on availability from BioMart web service.  Currently this argument can not be used}
\item{category}{Category is now replaced by page to better comply with the BioMart suite http://www.biomart.org}
\item{showGroups}{Availability of showGroups argument is pending on availability from BioMart web service.  Currently this argument can not be used}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
listAttributes(ensembl)
}
}
\keyword{methods}

