setMethod("show",signature(object="Mart"),
  function(object){	
    res = paste("Object of class 'Mart':\n Using the ",object@biomart," BioMart database\n Using the ",object@dataset," dataset\n", sep="")
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

setGeneric("martAttribPointers",def=function(obj,...)standardGeneric("martAttribPointers"))
setMethod("martAttribPointers",signature("Mart"),function(obj) obj@attributePointer)
setGeneric("martAttribPointers<-", function(obj, value) standardGeneric("martAttribPointers<-"))
setReplaceMethod("martAttribPointers","Mart",function(obj,value){
  obj@attributePointer <- value
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

setGeneric("martHost",def=function(obj,...)standardGeneric("martHost"))
setMethod("martHost",signature("Mart"), function(obj) obj@host)

setGeneric("martArchive",def=function(obj,...)standardGeneric("martArchive"))
setMethod("martArchive",signature("Mart"), function(obj) obj@archive)

setGeneric("martVSchema",def=function(obj,...)standardGeneric("martVSchema"))
setMethod("martVSchema",signature("Mart"), function(obj) obj@vschema)
setGeneric("martVSchema<-", function(obj, value) standardGeneric("martVSchema<-"))
setReplaceMethod("martVSchema","Mart",function(obj,value){
    obj@vschema <- value
      obj
  })



