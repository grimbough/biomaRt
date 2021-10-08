\name{listEnsemblArchives}
\alias{listEnsemblArchives}
\title{Lists the available archived versions of Ensembl}
\description{Returns a table containing the available archived versions of 
Ensembl, along with the dates they were created and the URL used to access
them.}

\usage{listEnsemblArchives(https)}

\arguments{
\item{https}{Deprecated argument.  Ensembl are enforcing https use from late 
2021 and this argument will be removed at this time as it no longer 
serves a purpose.  Originally - "Logical value of length 1.  Determines whether
https should be used to contact the Ensembl server."}
}

\author{Mike Smith} 

\examples{
listEnsemblArchives()
}
\keyword{methods}

