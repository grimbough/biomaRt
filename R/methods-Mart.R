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

setGeneric("martMainT",def=function(obj,...)standardGeneric("martMainT"))
setMethod("martMainT",signature("Mart"), function(obj) obj@mainTables)
setGeneric("martMainT<-", function(obj, value) standardGeneric("martMainT<-"))
setReplaceMethod("martMainT","Mart",function(obj,value){
  obj@mainTables <- value
  obj
})

setGeneric("martMySQL",def=function(obj,...)standardGeneric("martMySQL"))
setMethod("martMySQL",signature("Mart"),function(obj) obj@mysql)

setGeneric("martConnection",def=function(obj,...)standardGeneric("martConnection"))
setMethod("martConnection",signature("Mart"), function(obj) obj@connections)
setGeneric("martConnection<-", function(obj, value) standardGeneric("martConnection<-"))
setReplaceMethod("martConnection","Mart",function(obj,value){
  obj@connections <- value
  obj
})

setGeneric("martVSchema",def=function(obj,...)standardGeneric("martVSchema"))
setMethod("martVSchema",signature("Mart"), function(obj) obj@vschema)

setGeneric("martMySQLDriver",def=function(obj,...)standardGeneric("martMySQLDriver"))
setMethod("martMySQLDriver",signature("Mart"), function(obj) obj@mysqldriver)




