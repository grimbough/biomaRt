\name{NP2009code}
\alias{NP2009code}
\title{Display the analysis code from the 2009 Nature protocols paper}
\description{This function opens an editor displaying the analysis code of the Nature Protocols 2009 paper}
\usage{NP2009code()}

\details{The \code{\link{edit}} function uses \code{getOption("editor")} to select the editor.
Use, for instance, \code{options(editor="emacs")} to set another editor.}

\author{Steffen Durinck, Wolfgang Huber} 

\seealso{\code{\link{edit}}}

\examples{
if(interactive()){
NP2009code()
}
}
\keyword{methods}

