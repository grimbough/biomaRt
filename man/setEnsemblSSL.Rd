\name{setEnsemblSSL}
\alias{setEnsemblSSL}

\title{Save system specific SSL settings for contacting Ensembl}

\description{On some systems specific SSL settings have to be applied to allow
https connections to the Ensembl servers.  This function allows these to be 
saved in the biomaRt cache, so they will be retrieved each time they are 
needed.  biomaRt will try to determine them automatically, but this function
can be used to set them manually if required.}

\usage{
setEnsemblSSL(settings)
}

\arguments{
\item{settings}{A named list. Each entry should be a valid curl option,
as found in \code{\link[curl]{curl_options}}.}
}

\author{Mike Smith} 

\examples{

\dontrun{
  ssl_settings <- list("ssl_cipher_list" = "DEFAULT@SECLEVEL=1",
                       "ssl_verifypeer"  = FALSE)
  setEnsemblSSL
}
}


