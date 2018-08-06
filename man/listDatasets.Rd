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
    everything.}
}

\author{Steffen Durinck, Mike Smith} 

\examples{

if(interactive()){

#marts <- listMarts()
#index<-grep("ensembl",marts)

#mart <- useMart(marts[index])

#listDatasets(mart = mart)

#martDisconnect(mart = mart)
}
}
\keyword{methods}

