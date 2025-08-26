#' Class Mart
#'
#' Represents a Mart class, containing connections to different BioMarts
#'
#' @param object An object of class `Mart`
#'
#' @aliases Mart-class show,Mart-method
#' @section Methods: \describe{`show` Print summary of the object}
#' @author Steffen Durinck
#' @keywords methods
#' @rdname Mart-class
#' @importFrom methods setMethod signature
#' @export
setMethod("show", signature(object = "Mart"), function(object) {
  dbase <- ifelse(
    nzchar(object@biomart),
    yes = paste("  Using the", object@biomart, "BioMart database"),
    no = "  No database selected."
  )

  dset <- ifelse(
    nzchar(object@dataset),
    yes = paste("  Using the", object@dataset, "dataset"),
    no = "  No dataset selected."
  )

  res <- paste("Object of class 'Mart':", dbase, dset, sep = "\n")
  cat(res)
})

#' @importFrom methods setGeneric
setGeneric("martBM", function(obj, ...) standardGeneric("martBM"))
#' @importFrom methods setMethod signature
setMethod("martBM", signature("Mart"), function(obj) obj@biomart)
#' @importFrom methods setGeneric
setGeneric("martBM<-", function(obj, value) standardGeneric("martBM<-"))
#' @importFrom methods setReplaceMethod
setReplaceMethod("martBM", "Mart", function(obj, value) {
  obj@biomart <- value
  obj
})


#' @importFrom methods setGeneric
setGeneric("martAttributes", function(obj, ...) {
  standardGeneric("martAttributes")
})
#' @importFrom methods setMethod signature
setMethod("martAttributes", signature("Mart"), function(obj) obj@attributes)
#' @importFrom methods setGeneric
setGeneric("martAttributes<-", function(obj, value) {
  standardGeneric("martAttributes<-")
})
#' @importFrom methods setReplaceMethod
setReplaceMethod("martAttributes", "Mart", function(obj, value) {
  obj@attributes <- value
  obj
})


#' @importFrom methods setGeneric
setGeneric("martFilters", function(obj, ...) standardGeneric("martFilters"))
#' @importFrom methods setMethod signature
setMethod("martFilters", signature("Mart"), function(obj) obj@filters)
#' @importFrom methods setGeneric
setGeneric("martFilters<-", function(obj, value) {
  standardGeneric("martFilters<-")
})
#' @importFrom methods setReplaceMethod
setReplaceMethod("martFilters", "Mart", function(obj, value) {
  obj@filters <- value
  obj
})


#' @importFrom methods setGeneric
setGeneric("martDataset", function(obj, ...) standardGeneric("martDataset"))
#' @importFrom methods setMethod signature
setMethod("martDataset", signature("Mart"), function(obj) obj@dataset)
#' @importFrom methods setGeneric
setGeneric("martDataset<-", function(obj, value) {
  standardGeneric("martDataset<-")
})
#' @importFrom methods setReplaceMethod
setReplaceMethod("martDataset", "Mart", function(obj, value) {
  obj@dataset <- value
  obj
})

#' @importFrom methods setGeneric
setGeneric("martHost", function(obj, ...) standardGeneric("martHost"))
#' @importFrom methods setMethod signature
setMethod("martHost", signature("Mart"), function(obj) obj@host)
#' @importFrom methods setGeneric
setGeneric("martHost<-", function(obj, value) standardGeneric("martHost<-"))
#' @importFrom methods setReplaceMethod
setReplaceMethod("martHost", "Mart", function(obj, value) {
  obj@host <- value
  obj
})

#' @importFrom methods setGeneric
setGeneric("martVSchema", function(obj, ...) standardGeneric("martVSchema"))
#' @importFrom methods setMethod signature
setMethod("martVSchema", signature("Mart"), function(obj) obj@vschema)

#' @importFrom methods setGeneric
setGeneric("martHTTPConfig", function(obj) standardGeneric("martHTTPConfig"))
#' @importFrom methods setMethod signature
setMethod("martHTTPConfig", signature("Mart"), function(obj) {
  config <- do.call(c, obj@http_config)
  if (is.null(config)) {
    config <- list()
  }
  return(config)
})

#####################################################################
## new wrappers to enable keys, columns, select and keytypes

#' Retrieve information from the BioMart databases
#'
#' `select`, `columns` and `keys` are used together to extract
#' data from a `Mart` object.  These functions work much the same as the
#' classic biomaRt functions such as [getBM()] etc. and are provide here to
#' make this easier for people who are comfortable using these methods from
#' other Annotation packages.  Examples of other objects in other packages
#' where you can use these methods include (but are not limited to):
#' `ChipDb`, `OrgDb` `GODb`, `InparanoidDb` and
#' `ReactomeDb`.
#'
#' `columns` shows which kinds of data can be returned from the
#' `Mart` object.
#'
#' `keytypes` allows the user to discover which keytypes can be passed in
#' to `select` or `keys` as the `keytype` argument.
#'
#' `keys` returns keys from the `Mart` of the type specified by it's
#' `keytype` argument.
#'
#' `select` is meant to be used with these other methods and has arguments
#' that take the kinds of values that these other methods return.
#' `select` will retrieve the results as a data.frame based on parameters
#' for selected `keys` and `columns` and `keytype` arguments.
#'
#'
#' @aliases select-methods keys,Mart-method columns,Mart-method
#' keytypes,Mart-method select,Mart-method keys columns keytypes select
#' @param x the `Mart` object. The dataset of the `Mart` object must
#' already be specified for all of these methods.
#' @param keys the keys to select records for from the database.  Keys for some
#' keytypes can be extracted by using the `keys` method.
#' @param columns the columns or kinds of things that can be retrieved from the
#' database.  As with `keys`, all possible columns are returned by using
#' the `columns` method.
#' @param keytype the keytype that matches the keys used.  For the
#' `select` methods, this is used to indicate the kind of ID being used
#' with the keys argument. For the `keys` method this is used to indicate
#' which kind of keys are desired from `keys`
#' @param ... other arguments.  These include: \describe{ \item{pattern:}{the
#' pattern to match (used by keys)} \item{column:}{the column to search on.
#' This is used by keys and is for when the thing you want to pattern match is
#' different from the keytype, or when you want to simply want to get keys that
#' have a value for the thing specified by the column argument.}
#' \item{fuzzy:}{TRUE or FALSE value.  Use fuzzy matching? (this is used with
#' pattern by the keys method)} }
#' @return `keys`,`columns` and `keytypes` each return a
#' character vector or possible values.  `select` returns a data.frame.
#' @author Marc Carlson
#' @keywords methods
#' @examples
#'
#' if(interactive()) {
#'   ## 1st create a Mart object and specify the dataset
#'   mart <- useEnsembl(dataset="hsapiens_gene_ensembl",
#'                      biomart='ensembl')
#'   ## you can list the keytypes
#'   keytypes(mart)
#'   ## you can list the columns
#'   columns(mart)
#'   ## And you can extract keys when this is supported for your keytype of interest
#'   k = keys(mart, keytype="chromosome_name")
#'   head(k)
#'   ## You can even do some pattern matching on the keys
#'   k = keys(mart, keytype="chromosome_name", pattern="LRG")
#'   head(k)
#'   ## Finally you can use select to extract records for things that you are
#'   ## interested in.
#'   affy=c("202763_at","209310_s_at","207500_at")
#'   select(mart, keys=affy, columns=c('affy_hg_u133_plus_2','entrezgene_id'),
#'     keytype='affy_hg_u133_plus_2')
#'   }
#'
#' @importMethodsFrom AnnotationDbi keys keytypes columns select
#' @name select-methods

#' @rdname select-methods
#' @export
#' @importFrom methods setMethod
setMethod("keys", "Mart", function(x, keytype, ...) {
  # nolint next: undesirable_operator_linter.
  AnnotationDbi:::smartKeys(x = x, keytype = keytype, ..., FUN = .keys)
})
.keys <- function(x, keytype) {
  searchFilterOptions(mart = x, filter = keytype)
}

#' @rdname select-methods
#' @export
#' @importFrom methods setMethod
setMethod("keytypes", "Mart", function(x) listFilters(mart = x, what = "name"))
#' @rdname select-methods
#' @export
#' @importFrom methods setMethod
setMethod("columns", "Mart", function(x) {
  listAttributes(mart = x, what = "name")
})

## Arg checking is similar (but more limited) to what is done for getBM
#' @rdname select-methods
#' @export
#' @importFrom methods setMethod
setMethod("select", "Mart", function(x, keys, columns, keytype, ...) {
  if (missing(columns)) {
    stop("Argument 'columns' must be specified.")
  }
  if (!is.list(keytype) && keytype != "" && missing(keys)) {
    stop("Argument 'keys' must be specified.")
  }
  if (length(keytype) > 0 && length(keys) == 0) {
    stop("Keys argument contains no data.")
  }
  if (!(is.character(keytype)) || length(keytype) != 1) {
    stop("keytype should be single element character vector.")
  }
  getBM(attributes = columns, filters = keytype, values = keys, mart = x)
})
