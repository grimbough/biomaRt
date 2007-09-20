\name{listAttributes}
\alias{listAttributes}
\title{lists the attributes available in the selected dataset}
\description{Attributes are the outputs of a biomaRt query, they are the information we want to retrieve.  For example if we want to retrieve all entrez gene identifiers of genes located on chromosome X, entrezgene will be the attribute we use in the query.  The listAttributes function lists the available attributes in the selected dataset}
\usage{listAttributes(mart, group, category, showGroups = FALSE)}

\arguments{
\item{mart}{object of class Mart created using the useMart function}
\item{group}{Show only the attributes that belong to the specified attribute group.  To get a summary of the attribute groups set showGroups = TRUE or use the attributeSummary function}
\item{category}{Show only the attributes that belong to the specified attribute category.  To get a summary of the attribute pages set showGroups = TRUE or use the attributeSummary function}
\item{showGroups}{boolean to indicate if one wants to display the attribute categories and groups along with their names and descriptions}
}

\author{Steffen Durinck, http://www.stat.berkeley.edu/~steffen} 

\examples{

if(interactive()){
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
listAttributes(ensembl)
}
}
\keyword{methods}

