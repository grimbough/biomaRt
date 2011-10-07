\name{getXML}
\alias{getXML}
\title{Retrieves information from the BioMart database using an XML query}
\description{This function is a low level query function bypassing lots of biomaRts internal controls.  It allows for a direct XML query to a known BioMart webservice host.}
\usage{getXML(host="http://www.biomart.org/biomart/martservice?", xmlquery)}

\arguments{
  \item{host}{URL to BioMart webservice, is set to http://www.biomart.org/biomart/martservice? by default}
  \item{xmlquery}{XML query that needs to be send to the webservice}
}

\author{Steffen Durinck} 

\examples{
if(interactive()){
getXML(xmlquery="<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' uniqueRows = '1' count = '0' datasetConfigVersion = '0.6'> <Dataset name = 'hsapiens_gene_ensembl'><Attribute name = 'ensembl_gene_id'/><Filter name = 'chromosome_name' value = 'Y' /></Dataset></Query>")
}
}
\keyword{methods}

