setMethod("show",signature(object="Mart"),
          function(object){	
              
              dbase <- ifelse(nchar(object@biomart) != 0, 
                              yes = paste("  Using the", object@biomart, "BioMart database"),
                              no = "  No database selected.")
              
              dset <- ifelse(nchar(object@dataset) != 0, 
                             yes = paste("  Using the", object@dataset, "dataset"),
                             no = "  No dataset selected.")
              
              res <- paste("Object of class 'Mart':",
                           dbase,
                           dset,
                           sep="\n")
              cat(res)
          })

setGeneric("martBM",def=function(obj,...) standardGeneric("martBM"))
setMethod("martBM",signature("Mart"), function(obj) obj@biomart)
setGeneric("martBM<-", function(obj, value) standardGeneric("martBM<-"))
setReplaceMethod("martBM","Mart",function(obj,value){
  obj@biomart <- value
  obj
})


setGeneric("martAttributes",def=function(obj,...)standardGeneric("martAttributes"))
setMethod("martAttributes",signature("Mart"),function(obj) obj@attributes)
setGeneric("martAttributes<-", function(obj, value) standardGeneric("martAttributes<-"))
setReplaceMethod("martAttributes","Mart",function(obj,value){
  obj@attributes <- value
  obj
})


setGeneric("martFilters",def=function(obj,...)standardGeneric("martFilters"))
setMethod("martFilters",signature("Mart"),function(obj) obj@filters)
setGeneric("martFilters<-", function(obj, value) standardGeneric("martFilters<-"))
setReplaceMethod("martFilters","Mart",function(obj,value){
  obj@filters <- value
  obj
})


setGeneric("martDataset",def=function(obj,...)standardGeneric("martDataset"))
setMethod("martDataset",signature("Mart"), function(obj) obj@dataset)
setGeneric("martDataset<-", function(obj, value) standardGeneric("martDataset<-"))
setReplaceMethod("martDataset","Mart",function(obj,value){
  obj@dataset <- value
  obj
})

setGeneric("martHost", def=function(obj,...) standardGeneric("martHost"))
setMethod("martHost", signature("Mart"), function(obj) obj@host)
setGeneric("martHost<-", function(obj, value) standardGeneric("martHost<-"))
setReplaceMethod("martHost","Mart",function(obj,value){
  obj@host <- value
  obj
})

setGeneric("martVSchema",def=function(obj,...)standardGeneric("martVSchema"))
setMethod("martVSchema",signature("Mart"), function(obj) obj@vschema)

setGeneric("martHTTRConfig", def = function(obj) standardGeneric("martHTTRConfig"))
setMethod("martHTTRConfig", signature("Mart"), 
          function(obj) {
            config <- do.call(c, obj@httr_config)
            if(is.null(config)) {
              config <- httr::config()
            }
            return(config)
          }
)

#####################################################################
## new wrappers to enable keys, columns, select and keytypes
.keys <- function(x, keytype){
    searchFilterOptions(mart = x, filter = keytype)
}
setMethod("keys", "Mart",
    function(x, keytype, ...){
        AnnotationDbi:::smartKeys(x=x, keytype=keytype, ...,
                                  FUN=.keys)
    }
)

setMethod("keytypes", "Mart",
    function(x) listFilters(mart=x, what = "name")
)

setMethod("columns", "Mart",
    function(x) listAttributes(mart=x, what = "name")
)

## Arg checking is similar (but more limited) to what is done for getBM
setMethod("select", "Mart",
          function(x, keys, columns, keytype, ...){
              
              if(missing( columns ))
                  stop("Argument 'columns' must be specified.")              
              if(!is.list(keytype) && keytype != "" && missing( keys ))
                  stop("Argument 'keys' must be specified.")              
              if(length(keytype) > 0 && length(keys) == 0)
                  stop("Keys argument contains no data.")              
              if(!(is.character(keytype)) || length(keytype)!=1){
                  stop("keytype should be single element character vector.")
              }
              getBM(attributes=columns, filters=keytype, values=keys,  mart=x)
          }
)


