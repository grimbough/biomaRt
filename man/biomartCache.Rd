\name{biomartCache}
\alias{biomartCacheInfo}
\alias{biomartCacheClear}
\title{biomaRt result caching}
\description{biomaRt makes use of a results cache to speedup execution of
queries that have been run before. These functions provide details on the 
status of this cache, and allow it to be deleted.}

\usage{
biomartCacheInfo()
biomartCacheClear()
}

\value{These functions do not return anything and are called for their side 
effects.  \code{biomartCacheInfo()} prints the location of the cache, along
with the number of files and their total size on disk. 
\code{biomartCacheClear()} will delete the current contents of the cache.}

\author{Mike Smith} 

\keyword{IO}

