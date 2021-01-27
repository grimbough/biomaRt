\name{listEnsemblArchives}
\alias{listEnsemblArchives}
\title{Lists the available archived versions of Ensembl}
\description{Returns a table containing the available archived versions of 
Ensembl, along with the dates they were created and the URL used to access
them.}

\usage{listEnsemblArchives(https = FALSE)}

\arguments{
\item{https}{Logical value of length 1.  Determines whether https should be
used to contact the Ensembl server.  Defaults to \code{TRUE}.}
}

\author{Mike Smith} 

\examples{
listEnsemblArchives()
}
\keyword{methods}

