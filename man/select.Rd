%% actual page name
\alias{select-methods}
\name{select-methods}

%% for methods
\alias{keys,Mart-method}
\alias{columns,Mart-method}
\alias{keytypes,Mart-method}
\alias{select,Mart-method}

%% for convenience
\alias{keys}
\alias{columns}
\alias{keytypes}
\alias{select}


\title{Retrieve information from the BioMart databases}
\description{
  \code{select}, \code{columns} and \code{keys} are used together to
  extract data from a \code{Mart} object.  These functions work much the
  same as the classic biomaRt functions such as \code{getBM} etc. and
  are provide here to make this easier for people who are comfortable
  using these methods from other Annotation packages.  Examples of other
  objects in other packages where you can use these methods include (but
  are not limited to): \code{ChipDb}, \code{OrgDb} \code{GODb},
  \code{InparanoidDb} and \code{ReactomeDb}.

  \code{columns} shows which kinds of data can be returned from the
  \code{Mart} object.
  
  \code{keytypes} allows the user to discover which keytypes can be
  passed in to \code{select} or \code{keys} as the \code{keytype}
  argument.
  
  \code{keys} returns keys from the \code{Mart} of the type specified by
  it's \code{keytype} argument.
  
  \code{select} is meant to be used with these other methods and has
  arguments that take the kinds of values that these other methods
  return.  \code{select} will retrieve the results as a data.frame based
  on parameters for selected \code{keys} and \code{columns} and
  \code{keytype} arguments.
}
\usage{
  columns(x)
  keytypes(x)
  keys(x, keytype, ...)
  select(x, keys, columns, keytype, ...)
}

\arguments{
  \item{x}{the \code{Mart} object. The dataset of the \code{Mart} object
    must already be specified for all of these methods.}
  \item{keys}{the keys to select records for from the database.  Keys
    for some keytypes can be extracted by using the \code{keys} method.}
  \item{columns}{the columns or kinds of things that can be retrieved
    from the database.  As with \code{keys}, all possible columns are
    returned by using the \code{columns} method.}
  \item{keytype}{the keytype that matches the keys used.  For the
    \code{select} methods, this is used to indicate the kind of ID being used
    with the keys argument. For the \code{keys} method this is used to
    indicate which kind of keys are desired from \code{keys}
  }
  \item{...}{other arguments.  These include:
    \describe{
      \item{pattern:}{the pattern to match (used by keys)}
      \item{column:}{the column to search on. This is used by keys and is
        for when the thing you want to pattern match is different from
        the keytype, or when you want to simply want to get keys that
        have a value for the thing specified by the column argument.}
      \item{fuzzy:}{TRUE or FALSE value.  Use fuzzy matching? (this is
        used with pattern by the keys method)}
    }
  }
}

\value{
  \code{keys},\code{columns} and \code{keytypes} each return a character
    vector or possible values.  \code{select} returns a data.frame.
}

\author{Marc Carlson}

\examples{
if(interactive()) {
  ## 1st create a Mart object and specify the dataset
  mart <- useEnsembl(dataset="hsapiens_gene_ensembl",
                     biomart='ensembl')
  ## you can list the keytypes
  keytypes(mart)
  ## you can list the columns
  columns(mart)
  ## And you can extract keys when this is supported for your keytype of interest
  k = keys(mart, keytype="chromosome_name")
  head(k)
  ## You can even do some pattern matching on the keys
  k = keys(mart, keytype="chromosome_name", pattern="LRG")
  head(k)
  ## Finally you can use select to extract records for things that you are
  ## interested in.
  affy=c("202763_at","209310_s_at","207500_at")
  select(mart, keys=affy, columns=c('affy_hg_u133_plus_2','entrezgene_id'),
    keytype='affy_hg_u133_plus_2')
  }
}
\keyword{methods}
