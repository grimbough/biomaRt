\name{cleanBM}
\alias{cleanBM}
\title{Clear getBM result from duplicates}
\description{This function removes duplicated results retrieved by
  getBM.  Two sets of IDs are required in the query result to be able to
  use this function.  The query id that was used in  the getBM function
  as filter and was also used as attribute to retrieve
  e.g. affy\_hg\_u95av2 and a result id e.g. GO identifiers that were retrieved.}
\usage{cleanBM(bmresult, query_id_col=1, result_id_col=2)}

\arguments{
  \item{bmresult}{A data.frame that was returned by the getBM function}
  \item{query_id_col}{Index of column of the query id in the result}
  \item{result_id_col}{Index of column of the result id in the result}
}

\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\examples{
if(interactive()){
mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart<-useDataset("hsapiens_gene_ensembl",mart)
bm=getBM(attributes=c("affy_hg_u95av2","go","hgnc_symbol","chromosome_name","band"),filters="affy_hg_u95av2",values=c("1939_at","1000_at","1454_at"), mart=mart)

writeLines("## before duplicate removal ##")

print(bm)

writeLines("## After duplicate removal ##")
print(cleanBM(bm,1,2))
}
}
\keyword{methods}

