\name{getXML}
\alias{getXML}
\title{Retrieves information from the BioMart database using an XML query}
\description{DEPRECATED! - This function is a low level query function bypassing lots of biomaRts internal controls.  It allows for a direct XML query to a known BioMart webservice host.}

\usage{getXML(host="http://www.ensembl.org/biomart/martservice?", xmlquery)}

\arguments{
  \item{host}{URL to BioMart webservice, is set to http://www.ensembl.org/biomart/martservice? by default}
  \item{xmlquery}{XML query that needs to be send to the webservice}
}

\details{This function has now been deprecated.  Similar functionality can be 
achieved using the unexported function `\code{biomaRt:::.submitQueryXML}.  If 
you are still using this function and this affects your work, please contact 
the package maintainer at https://github.com/grimbough/biomaRt/issues}

\author{Steffen Durinck, Mike Smith} 

\examples{
if(interactive()){
xmlquery="<?xml version='1.0' encoding='UTF-8'?>
            <!DOCTYPE Query>
            <Query virtualSchemaName = 'default' uniqueRows = '1' 
                   count = '0' datasetConfigVersion = '0.6'>
                <Dataset name = 'hsapiens_gene_ensembl'>
                    <Attribute name = 'ensembl_gene_id'/>
                    <Filter name = 'chromosome_name' value = 'Y' />
                </Dataset>
            </Query>"
getXML(host = "www.ensembl.org/biomart/martservice?", xmlquery = xmlquery)
}
}
\keyword{methods}

