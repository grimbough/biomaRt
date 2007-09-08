setMethod("show",signature(object="Mart"),
  function(object){	
    res = paste("Object of class 'Mart':\n Using the ",object@biomart," BioMart database\n Using the ",object@dataset," dataset\n", sep="")
    cat(res)
})

setGeneric("martBM",def=function(obj,...) standardGeneric("martBM"))
setMethod("martBM",signature("Mart"), function(obj) obj@biomart)

setGeneric("martAttributes",def=function(obj,...)standardGeneric("martAttributes"))
setMethod("martAttributes",signature("Mart"),function(obj) obj@attributes)

setGeneric("martAttribPointers",def=function(obj,...)standardGeneric("martAttribPointers"))
setMethod("martAttribPointers",signature("Mart"),function(obj) obj@attributePointer)

setGeneric("martFilters",def=function(obj,...)standardGeneric("martFilters"))
setMethod("martFilters",signature("Mart"),function(obj) obj@filters)

setGeneric("martDataset",def=function(obj,...)standardGeneric("martDataset"))
setMethod("martDataset",signature("Mart"), function(obj) obj@dataset)

setGeneric("martHost",def=function(obj,...)standardGeneric("martHost"))
setMethod("martHost",signature("Mart"), function(obj) obj@host)

setGeneric("martMainT",def=function(obj,...)standardGeneric("martMainT"))
setMethod("martMainT",signature("Mart"), function(obj) obj@mainTables)

setGeneric("martMySQL",def=function(obj,...)standardGeneric("martMySQL"))
setMethod("martMySQL",signature("Mart"),function(obj) obj@mysql)

setGeneric("martConnection",def=function(obj,...)standardGeneric("martConnection"))
setMethod("martConnection",signature("Mart"), function(obj) obj@connections)

setGeneric("martVSchema",def=function(obj,...)standardGeneric("martVSchema"))
setMethod("martVSchema",signature("Mart"), function(obj) obj@vschema)





