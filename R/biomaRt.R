##########################
#biomaRt source code     #
##########################
#                        #
#Licence: Artistic       #
#Author: Steffen Durinck #
##########################

###############
#messageToUser#
###############

messageToUser <- function(message){
 if(interactive()){
  cat(message)
 }
}

##############################################################
#martCheck                                                   # 
#                                                            #
#This function checks if there is a valid Mart object,       # 
#if a dataset is selected and                                #
#if the correct BioMart database has been selected (optional)# 
##############################################################

martCheck = function(mart, biomart = NULL){
  if( missing( mart ) || class( mart ) != 'Mart'){
    stop("You must provide a valid Mart object. To create a Mart object use the function: useMart.  Check ?useMart for more information.")
  }
  if(!is.null(biomart)){
    martcheck = strsplit(martBM(mart),"_", fixed = TRUE, useBytes = TRUE)[[1]][1]
    if(martcheck[1] != biomart)stop(paste("This function only works when used with the ",biomart," BioMart.",sep="")) 
  }
  if(martDataset(mart)==""){
    stop("No dataset selected, please select a dataset first.  You can see the available datasets by using the listDatasets function see ?listDatasets for more information.  Then you should create the Mart object by using the useMart function.  See ?useMart for more information");
  }
}

checkWrapperArgs = function(id, type, mart){
 if(missing(type))stop("Specify the type of identifier you are using, see ?getGene for details. Valid values for the type argument can be found with the listFilters function.") 
  if(!type %in% listFilters(mart)[,1])stop(paste("Invalid identifier type:",type," see ?getGene for details. Use the listFilters function to get the valid value for the type argument.", sep=""))   
  if(missing(id))stop("No identifiers specified.  Use the id argument to specify a vector of identifiers for which you want to retrieve the annotation.")
}
#######################################################
#listMarts:                                           #
#list all available BioMart databases by default      #
#listMarts will check the central service to see which#
#BioMart databases are present                        #
#######################################################

bmRequest <- function(request, ssl.verifypeer = TRUE, verbose = FALSE){
  if(verbose) writeLines(paste("Attempting web service request:\n",request, sep=""))
  result = tryCatch(getURL(request, ssl.verifypeer = ssl.verifypeer,followlocation = TRUE), error = function(e){ cat("Request to BioMart web service failed. Verify if you are still connected to the internet.  Alternatively the BioMart web service is temporarily down.  Check http://www.biomart.org and verify if this website is available.\n")})
  return(result)
}

listMarts <- function( mart, host="www.biomart.org", path="/biomart/martservice", port=80,includeHosts = FALSE, archive = FALSE, ssl.verifypeer = TRUE, verbose = FALSE){

  if(archive){
    request = paste("http://",host,":",port,path,"?type=registry_archive&requestid=biomaRt", sep="")
    registry = bmRequest(request = request, ssl.verifypeer = ssl.verifypeer, verbose = verbose)
  }
  else{
    request = paste("http://",host,":",port,path,"?type=registry&requestid=biomaRt", sep="")	
    registry = bmRequest(request = request, ssl.verifypeer = ssl.verifypeer, verbose = verbose)
  }
  registry = xmlTreeParse(registry, asText=TRUE)
  registry = registry$doc$children[[1]]
  
  marts = list(biomart = NULL, version = NULL, host = NULL, path = NULL, database = NULL)
  index = 1
  
  if(host != "www.biomart.org" || archive){
    for(i in seq(len=xmlSize(registry))){
      if(xmlName(registry[[i]])=="MartURLLocation"){  
        if(xmlGetAttr(registry[[i]],"visible") == 1){
          if(!is.null(xmlGetAttr(registry[[i]],"name"))) marts$biomart[index] = xmlGetAttr(registry[[i]],"name")
          if(!is.null(xmlGetAttr(registry[[i]],"database"))) marts$database[index] = xmlGetAttr(registry[[i]],"database")
          if(!is.null(xmlGetAttr(registry[[i]],"displayName"))) marts$version[index] = xmlGetAttr(registry[[i]],"displayName")
          if(!is.null(xmlGetAttr(registry[[i]],"host"))) marts$host[index] = xmlGetAttr(registry[[i]],"host")
          if(!is.null(xmlGetAttr(registry[[i]],"path"))) marts$path[index] = xmlGetAttr(registry[[i]],"path")
          if(!is.null(xmlGetAttr(registry[[i]],"port"))) marts$port[index] = xmlGetAttr(registry[[i]],"port")
          if(!is.null(xmlGetAttr(registry[[i]],"serverVirtualSchema"))){
            marts$vschema[index] =  xmlGetAttr(registry[[i]],"serverVirtualSchema")
          }
          index=index+1
        }
      }
    }
  }
  else{
    for(i in seq(len=xmlSize(registry))){
      if(xmlName(registry[[i]])=="MartURLLocation"){  
        if(xmlGetAttr(registry[[i]],"visible") == 1){
          if(!is.null(xmlGetAttr(registry[[i]],"name"))) marts$biomart[index] = xmlGetAttr(registry[[i]],"name")
          if(!is.null(xmlGetAttr(registry[[i]],"database"))) marts$database[index] = xmlGetAttr(registry[[i]],"database")
          if(!is.null(xmlGetAttr(registry[[i]],"displayName"))) marts$version[index] = xmlGetAttr(registry[[i]],"displayName")
          marts$host[index] = host
          marts$path[index] = path
          marts$port[index] = 80
          if(!is.null(xmlGetAttr(registry[[i]],"serverVirtualSchema"))){
            marts$vschema[index] =  xmlGetAttr(registry[[i]],"serverVirtualSchema")
          }
          index=index+1
        }
      }
    }
  }
  if(includeHosts){
    return(marts)
  }
  else{
    if(archive){
      ret = data.frame(marts$database,marts$version)
    }
    else{
      ret = data.frame(marts$biomart,marts$version)
    }
    colnames(ret)=c("biomart","version")
    return(ret)
  } 
}

#################################
# #                           # #
# # Generic BioMart functions # #
# #                           # #
#################################

useMart <- function(biomart, dataset, host = "www.biomart.org", path = "/biomart/martservice", port = 80, archive = FALSE, ssl.verifypeer = TRUE, version, verbose = FALSE){

  if(missing(biomart) && missing(version)) stop("No biomart databases specified. Specify a biomart database to use using the biomart or version argument")
  if(!missing(biomart)){ 
  if(!(is.character(biomart)))
      stop("biomart argument is no string.  The biomart argument should be a single character string")
  }

  marts=NULL
  marts=listMarts(host=host, path=path, port=port, includeHosts = TRUE, archive = archive, ssl.verifypeer = ssl.verifypeer)
  mindex = NA
  if(!missing(biomart)){ 
   mindex=match(biomart,marts$biomart)
  }
  if(!missing(version)){
   mindex=match(version,marts$version)
  }
  if(is.na(mindex) || archive){
    mindex=match(biomart,marts$database)
  }
  if(is.na(mindex))
    stop("Incorrect BioMart name, use the listMarts function to see which BioMart databases are available")
    
    if(is.na(marts$path[mindex]) || is.na(marts$vschema[mindex]) || is.na(marts$host[mindex]) || is.na(marts$port[mindex]) || is.na(marts$path[mindex])) stop("The selected biomart databases is not available due to error in the BioMart central registry, please report so the BioMart registry file can be fixed.")
    if(marts$path[mindex]=="") marts$path[mindex]="/biomart/martservice" #temporary to catch bugs in registry
   if(archive) biomart = marts$biomart[mindex]
   if(!missing(version)) biomart = marts$biomart[mindex]
    biomart = sub(" ","%20",biomart, fixed = TRUE, useBytes = TRUE)
    mart <- new("Mart", biomart = biomart,vschema = marts$vschema[mindex], host = paste("http://",marts$host[mindex],":",marts$port[mindex],marts$path[mindex],sep=""), archive = archive)
    BioMartVersion=bmVersion(mart, verbose=verbose)
    if(martHost(mart) =="http://www.biomart.org:80/biomart/martservice"){
      if(verbose) writeLines("Using Central Repository at www.biomart.org");
      martVSchema(mart) <- 'default'  #Assume central service query uses default vSchema 
    }
    if(verbose){
      writeLines(paste("BioMartServer running BioMart version:",BioMartVersion,sep=" "))
      writeLines(paste("Mart virtual schema:",martVSchema(mart),sep=" "))
      writeLines(paste("Mart host:",martHost(mart),sep=" "))
    }
    if(!missing(dataset)){
      mart = useDataset(mart = mart, dataset=dataset, verbose = verbose)
    }
    return(mart)
}

listDatasets <- function(mart, verbose = FALSE) {
  if(missing(mart) || !is(mart, 'Mart'))
    stop("No Mart object given or object not of class 'Mart'")
  request = paste(martHost(mart),"?type=datasets&requestid=biomaRt&mart=",martBM(mart),sep="")
  bmResult = bmRequest(request = request, verbose = verbose)
  con = textConnection(bmResult)
  txt = scan(con, sep="\t", blank.lines.skip=TRUE, what="character", quiet=TRUE)
  #txt = tryCatch(scan(request, sep="\t", blank.lines.skip=TRUE, what="character", quiet=TRUE), error = function(e){stop("Request to BioMart web service failed. Verify if you are still connected to the internet.  Alternatively the BioMart web service is temporarily down.")})
  close(con)
   
  ## select visible ("1") table sets
    i = intersect(which(txt=="TableSet"), which(txt=="1")-3L)
  
  res = data.frame(dataset     = I(txt[i+1L]),
    description = I(txt[i+2L]),
    version     = I(txt[i+4L]))
  return(res)
}

## Check version of BioMart service

bmVersion <- function(mart, verbose=FALSE){
  request = ""
  request = paste(martHost(mart),"?type=version","&requestid=biomaRt&mart=",martBM(mart),sep="")
  BioMartVersion = bmRequest(request = request, verbose = verbose)
  con = textConnection(BioMartVersion)
  bmVersionParsed = read.table(con, sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
  close(con)
  if(verbose) print(bmVersionParsed)
  bmv=""
  if(dim(bmVersionParsed)[2] >=1){
    bmv=bmVersionParsed[1,1]
  }
  return(bmv)
}

## Retrieve attributes and filters from web service

bmAttrFilt <- function(type, mart, verbose=FALSE){
  request = ""
  request = paste(martHost(mart),"?type=",type,"&dataset=",martDataset(mart),"&requestid=biomaRt&mart=",martBM(mart),"&virtualSchema=",martVSchema(mart),sep="")
  attrfilt = bmRequest(request = request, verbose = verbose)
  con = textConnection(attrfilt)
  attrfiltParsed = read.table(con, sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
  close(con)
  if(type=="attributes"){
    if(dim(attrfiltParsed)[2] < 3)
      stop("biomaRt error: looks like we're connecting to incompatible version of BioMart suite.")
    cnames = seq_len(dim(attrfiltParsed)[2])
    cnames=paste(type,cnames,sep="")
    cnames[1] = "name"
    cnames[2] = "description"
    cnames[3] = "fullDescription"
    if(dim(attrfiltParsed)[2] < 4){
      warning("biomaRt warning: looks like we're connecting to an older version of BioMart suite. Some biomaRt functions might not work.")
    }
    else{
      cnames[4] = "page"
    }
    colnames(attrfiltParsed) = cnames
  }
   if(type=="filters"){
     if(dim(attrfiltParsed)[2] < 4)
       stop("biomaRt error: looks like we're connecting to incompatible version of BioMart suite.")
     cnames = seq(1:dim(attrfiltParsed)[2])
     cnames=paste(type,cnames,sep="")
     cnames[1] = "name"
     cnames[2] = "description"
     cnames[3] = "options"
     cnames[4] = "fullDescription"
     if(dim(attrfiltParsed)[2] < 7){
       warning("biomaRt warning: looks like we're connecting to an older version of BioMart suite. Some biomaRt functions might not work.")
     }
     else{
       cnames[5] = "filters"
       cnames[6] = "type"
       cnames[7] = "operation"
     }
     colnames(attrfiltParsed) = cnames
   }
  return(attrfiltParsed)
}

## Select a BioMart dataset             

useDataset <- function(dataset, mart, verbose = FALSE){
  if(missing(mart) || class(mart)!="Mart") stop("No valid Mart object given, specify a Mart object with the attribute mart")
  if(missing(dataset)) stop("No dataset given.  Please use the dataset argument to specify which dataset you want to use. Correct dataset names can be obtained with the listDatasets function.")
validDatasets=listDatasets(mart)
  if(is.na(match(dataset, validDatasets$dataset)))stop(paste("The given dataset: ",dataset,", is not valid.  Correct dataset names can be obtained with the listDatasets function."))
  martDataset(mart) = dataset  
  if(verbose) messageToUser("Checking attributes ...")
  martAttributes(mart) <- bmAttrFilt("attributes",mart, verbose = verbose)
  if(verbose){
    messageToUser(" ok\n")
    messageToUser("Checking filters ...")
  }
  martFilters(mart) <- bmAttrFilt("filters",mart, verbose = verbose)
  if(verbose) messageToUser(" ok\n")
  return( mart )
}

## getName

getName <- function(x, pos) if(is.null(x[[pos]])) NA else x[[pos]] 

## listAttributes

listAttributes <- function(mart, page, what = c("name","description"), group, category, showGroups = FALSE) {
  martCheck(mart)
  if(!missing(group)){
    .Defunct(msg="Currently the group argument is defunct.  Pending on availability from the BioMart web service group will become activated again or not.")
  }
  if(showGroups){
    .Defunct(msg="Currently the showGroups argument is defunct.  Pending on availability from the BioMart web service showGroups will become activated again or not.")
  }
  
  if(!missing(category)){
    page = category
    warning("To better comply with the BioMart suite (http://www.biomart.org) attribute categories are now known as attribute pages and can be specified with the page parameter")
  }
  if(!missing(page) && !page %in% attributePages(mart)) stop(paste("The chosen page: ",page," is not valid, please use the correct page name using the attributePages function",sep=""))
  attrib=NULL
  if(!missing(page)){
    sel = which(martAttributes(mart)[,"page"] == page)
    attrib = martAttributes(mart)[sel,what]
  }
  else{
    attrib = martAttributes(mart)[,what]
  }
  return(attrib)
}

## attributePages

attributePages <- function(mart){
  martCheck(mart)
  pages = unique(martAttributes(mart)[,"page"])
  return(pages)
}

## listFilters

listFilters <- function(mart, what = c("name", "description"), group = "DEFUNCT") {

  if(!missing(group))
    .Defunct(msg = "The argument 'group' is defunct. If you need advice how to replace that functionality, please contact the package maintainer for advice.")
  
  martCheck(mart)
  filters = martFilters(mart)
  badwhat = !(what %in% colnames(filters))
  if(any(badwhat))
    stop(sprintf("The function argument 'what' contains %s: %s\nValid are: %s\n",
                 if(sum(badwhat)>1) "invalid values" else "an invalid value",
                 paste(what[badwhat], collapse=", "),
                 paste(colnames(filters), collapse=", ")))
  return(filters[, what])
}

## filterOptions

filterOptions <- function(filter, mart){
  if(missing(filter)) stop("No filter given. Please specify the filter for which you want to retrieve the possible values.")
  if(class(filter)!="character")stop("Filter argument should be of class character")
  martCheck(mart)
  if(!filter %in% listFilters(mart, what="name"))stop("Filter not valid, check for typo in filter argument.")
  sel = which(listFilters(mart, what="name") == filter)
  return(listFilters(mart,what="options")[sel])
}

## filterType

filterType <- function(filter, mart){
  if(missing(filter)) stop("No filter given. Please specify the filter for which you want to retrieve the filter type")
  if(class(filter)!="character")stop("Filter argument should be of class character")
  martCheck(mart)
  type="unknown"
  sel = which(listFilters(mart, what="name") == filter)
  if(is.null(sel))stop(paste("Invalid filter",filter, sep=": "))
  type = listFilters(mart,what="type")[sel]
  return(type)
}

##########################################
#getBM: generic BioMart query function   # 
##########################################

getBM <- function(attributes, filters = "", values = "", mart, curl = NULL, checkFilters = TRUE, verbose=FALSE, uniqueRows=TRUE, bmHeader=TRUE){
  
  martCheck(mart)
  if(missing( attributes ))
    stop("Argument 'attributes' must be specified.")
  
   if(is.list(filters) && !missing( values ))
             warning("Argument 'values' should not be used when argument 'filters' is a list and will be ignored.")
   if(is.list(filters) && is.null(names(filters)))
             stop("Argument 'filters' must be a named list when sent as a list.")
   if(!is.list(filters) && filters != "" && missing( values ))
        stop("Argument 'values' must be specified.")

  if(length(filters) > 0 && length(values) == 0)
    stop("Values argument contains no data.")
  
  if(is.list(filters)){
    values = filters
    filters = names(filters)
  }
  
  if(class(uniqueRows) != "logical")
    stop("Argument 'uniqueRows' must be a logical value, so either TRUE or FALSE")

  xmlQuery = paste("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = '",martVSchema(mart),"' uniqueRows = '",as.numeric(uniqueRows),"' count = '0' datasetConfigVersion = '0.6' header='",as.numeric(bmHeader),"' requestid= \"biomaRt\"> <Dataset name = '",martDataset(mart),"'>",sep="")
  
  #checking the Attributes
  invalid = !(attributes %in% listAttributes(mart, what="name"))
  if(any(invalid))
    stop(paste("Invalid attribute(s):", paste(attributes[invalid], collapse=", "),
               "\nPlease use the function 'listAttributes' to get valid attribute names"))

  #check if attributes come from multiple attribute pages currently disabled until ID issue resovled at Ensembl
  if(FALSE){
    att = listAttributes(mart, what=c("name","page"))
    att = att[which(att[,1] %in% attributes),]
    attOK = FALSE
    pages = unique(att[,2])
    if(length(pages) <= 1){
      attOK = TRUE
    }
    else{
      for(page in pages){
        if(length(attributes) == length(which(attributes %in% att[which(att[,2] == page),1]))) attOK = TRUE
      }
    }
    if(!attOK){
      stop(paste("Querying attributes from multiple attribute pages is not allowed.  To see the attribute pages attributes belong to, use the function attributePages."))
    }
  }
  #attribute are ok lets add them to the query
  attributeXML =  paste("<Attribute name = '", attributes, "'/>", collapse="", sep="")
  
  #checking the filters
  if(filters[1] != "" && checkFilters){
    invalid = !(filters %in% listFilters(mart, what="name"))
    if(any(invalid))
      stop(paste("Invalid filters(s):", paste(filters[invalid], collapse=", "),
                 "\nPlease use the function 'listFilters' to get valid filter names"))
  }
  
  filterXML = NULL
  
  if(length(filters) > 1){
    if(class(values)!= "list")stop("If using multiple filters, the 'value' has to be a list.\nFor example, a valid list for 'value' could be: list(affyid=c('1939_at','1000_at'), chromosome= '16')\nHere we select on Affymetrix identifier and chromosome, only results that pass both filters will be returned");
   
    for(i in seq(along = filters)){
      if(filters[i] %in% listFilters(mart, what = "name")){
        filtertype=filterType(filters[i], mart)
        if(filtertype == 'boolean' || filtertype == 'boolean_list'){
          if(!is.logical(values[[i]])) stop(paste("biomaRt error: ",filters[i]," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="")  
          if(!values[[i]]){
            values[[i]] = 1
          }
          else{
            values[[i]] = 0 
          }
          filterXML = paste(filterXML,paste("<Filter name = '",filters[i],"' excluded = \"",values[[i]],"\" />", collapse="",sep=""),sep="")
        }
        else{
          if(is.numeric(values[[i]])){ values[[i]] = as.integer(values[[i]])}
          valuesString = paste(values[[i]],"",collapse=",",sep="")
          filterXML = paste(filterXML,paste("<Filter name = '",filters[i],"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
        }
      }
      else{ #used for attributes with values as these are treated as filters in BioMart
        valuesString = paste(values[[i]],"",collapse=",",sep="")
        filterXML = paste(filterXML,paste("<Filter name = '",filters[i],"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
      } 
    }
  }
  else{
    if(filters != ""){
      if(is.list(values)){
        values = unlist(values)
      }
      if(filters %in% listFilters(mart, what="name")){
        filtertype =filterType(filters, mart)
        if(filtertype == 'boolean' || filtertype == 'boolean_list'){
          if(!is.logical(values)) stop(paste("biomaRt error: ",filters," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="") 
          if(!values){
            values = 1
          }
          else{
            values = 0 
          }
          filterXML = paste("<Filter name = '",filters,"' excluded = \"",values,"\" />", collapse="",sep="")
        }
        else{
          if(is.numeric(values)){
            values = as.integer(values)
          }  
          valuesString = paste(values,"",collapse=",",sep="")
          filterXML = paste("<Filter name = '",filters,"' value = '",valuesString,"' />", collapse="",sep="")
        }
      }
      else{ #used for attributes with values as these are treated as filters in BioMart
        valuesString = paste(values,"",collapse=",",sep="")
        filterXML = paste(filterXML,paste("<Filter name = '",filters,"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
      }
    }
    else{
      filterXML=""
    }
  }
  
  xmlQuery = paste(xmlQuery, attributeXML, filterXML,"</Dataset></Query>",sep="")
    
  if(verbose){
    cat(paste(xmlQuery,"\n", sep=""))
  }      

  postRes = tryCatch(postForm(paste(martHost(mart),"?",sep=""),"query" = xmlQuery), error = function(e){stop("Request to BioMart web service failed. Verify if you are still connected to the internet.  Alternatively the BioMart web service is temporarily down.")})
  
  if(!(is.character(postRes) && (length(postRes)==1L)))
    stop("The query to the BioMart webservice returned an invalid result: biomaRt expected a character string of length 1. Please report this to the mailing list.")

  if(gsub("\n", "", postRes, fixed = TRUE, useBytes = TRUE) == "") {
    
    result = as.data.frame(matrix("", ncol=length(attributes), nrow=0), stringsAsFactors=FALSE)
    
  } else {
    
    if(length(grep("^Query ERROR", postRes))>0L)
      stop(postRes)

    ## convert the serialized table into a dataframe
    con = textConnection(postRes)
    result = read.table(con, sep="\t", header=bmHeader, quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors=FALSE)
    close(con)

    if(!(is(result, "data.frame") && (ncol(result)==length(attributes)))) {
      print(head(result))
      stop("The query to the BioMart webservice returned an invalid result: the number of columns in the result table does not equal the number of attributes in the query. Please report this to the mailing list.")
    }
  }
  if(!bmHeader){
    colnames(result) = attributes
  }
  else{
    att = listAttributes(mart)
    colm = match(colnames(result),att[,2])
    if(length(which(is.na(colm))) > 0){
      warning("Unable to match column names of BioMart output")
    }
    else{
      colnames(result) = att[colm,1]
      result = result[,attributes]
    }
  }
  return(result)
}

###################################
#getLDS: Multiple dataset linking #
###################################

getLDS <- function(attributes, filters = "", values = "", mart, attributesL, filtersL = "", valuesL = "", martL, verbose = FALSE, uniqueRows = TRUE, bmHeader = TRUE) {
  
  martCheck(mart)
  martCheck(martL)
  invalid = !(attributes %in% listAttributes(mart, what="name"))
  if(any(invalid))
    stop(paste("Invalid attribute(s):", paste(attributes[invalid], collapse=", "),
               "\nPlease use the function 'listAttributes' to get valid attribute names"))
  
  invalid = !(attributesL %in% listAttributes(martL, what="name"))
  if(any(invalid))
    stop(paste("Invalid attribute(s):", paste(attributesL[invalid], collapse=", "),
               "\nPlease use the function 'listAttributes' to get valid attribute names"))
  
  if(filters[1] != ""){
    invalid = !(filters %in% listFilters(mart, what="name"))
    if(any(invalid))
      stop(paste("Invalid filters(s):", paste(filters[invalid], collapse=", "),
                 "\nPlease use the function 'listFilters' to get valid filter names"))
  }
  if(filtersL[1] != ""){
    invalid = !(filtersL %in% listFilters(martL, what="name"))
    if(any(invalid))
      stop(paste("Invalid filters(s):", paste(filtersL[invalid], collapse=", "),
                 "\nPlease use the function 'listFilters' to get valid filter names"))
  }
  
  xmlQuery = paste("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' uniqueRows = '",as.numeric(uniqueRows),"' count = '0' datasetConfigVersion = '0.6' header='",as.numeric(bmHeader),"' requestid= \"biomaRt\"> <Dataset name = '",martDataset(mart),"'>",sep="")
  attributeXML = paste("<Attribute name = '", attributes, "'/>", collapse="", sep="")
  if(length(filters) > 1){
    if(class(values)!= "list")
      stop("If using multiple filters, the 'value' has to be a list.\nFor example, a valid list for 'value' could be: list(affyid=c('1939_at','1000_at'), chromosome= '16')\nHere we select on affyid and chromosome, only results that pass both filters will be returned");
    filterXML = NULL
    for(i in seq(along=filters)){
      if(filterType(filters[i],mart) == 'boolean' || filterType(filters[i],mart) == 'boolean_list'){
        if(!is.logical(values[[i]])) stop(paste("biomaRt error: ",filters[i]," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="") 
        if(!values[[i]]){
          values[[i]] = 1
        }
        else{
          values[[i]] = 0 
        }
        filterXML = paste(filterXML,paste("<Filter name = '",filters[i],"' excluded = \"",values[[i]],"\" />", collapse="",sep=""),sep="")
      }
      else{
        valuesString = paste(values[[i]],"",collapse=",",sep="")
        filterXML = paste(filterXML,paste("<Filter name = '",filters[i],"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
      }
    }
  }
  else{
    if(filters != ""){       
      if(filterType(filters,mart) == 'boolean' || filterType(filters,mart) == 'boolean_list'){
        if(!is.logical(values)) stop(paste("biomaRt error: ",filters," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="") 
        if(!values){
          values = 1
        }
        else{
          values = 0 
        }
        
        filterXML = paste("<Filter name = '",filters,"' excluded = \"",values,"\" />", collapse="",sep="")
      }
      else{
        valuesString = paste(values,"",collapse=",",sep="")
        filterXML = paste("<Filter name = '",filters,"' value = '",valuesString,"' />", collapse="",sep="")
      }
    }
    else{
      filterXML=""
    }
  }
  
  xmlQuery = paste(xmlQuery, attributeXML, filterXML,"</Dataset>",sep="")
  xmlQuery = paste(xmlQuery, "<Dataset name = '",martDataset(martL),"' >", sep="")
  linkedAttributeXML =  paste("<Attribute name = '", attributesL, "'/>", collapse="", sep="")
  
  if(length(filtersL) > 1){
    if(class(valuesL)!= "list")
      stop("If using multiple filters, the 'value' has to be a list.\nFor example, a valid list for 'value' could be: list(affyid=c('1939_at','1000_at'), chromosome= '16')\nHere we select on affyid and chromosome, only results that pass both filters will be returned");
    linkedFilterXML = NULL
    for(i in seq(along=filtersL)){
      if(filterType(filtersL,martL) == 'boolean' || filterType(filtersL,martL) == 'boolean_list'){
        if(!is.logical(valuesL[[i]])) stop(paste("biomaRt error: ",filtersL[i]," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="") 
        if(!valuesL[[i]]){
          valuesL[[i]] = 1
        }
        else{
          valuesL[[i]] = 0 
        } 
        linkedFilterXML = paste(linkedFilterXML,paste("<Filter name = '",filtersL[i],"' excluded = \"",valuesL[[i]],"\" />", collapse="",sep=""),sep="")
      }
      else{
        valuesString = paste(valuesL[[i]],"",collapse=",",sep="")
        linkedFilterXML = paste(linkedFilterXML,paste("<Filter name = '",filtersL[i],"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
      }
    }
  }
  else{
    if(filtersL != ""){
      if(filterType(filtersL,martL) == 'boolean' || filterType(filtersL,martL) == 'boolean_list'){
        if(!is.logical(valuesL)) stop(paste("biomaRt error: ",filtersL," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="") 
        if(!valuesL){
          valuesL = 1
        }
        else{
          valuesL = 0 
        }
        
        linkedFilterXML = paste("<Filter name = '",filtersL,"' excluded = \"",valuesL,"\" />", collapse="",sep="")
      }
      else{
        valuesString = paste(valuesL,"",collapse=",",sep="")
        linkedFilterXML = paste("<Filter name = '",filtersL,"' value = '",valuesString,"' />", collapse="",sep="")
      }
    }
    else{
      linkedFilterXML=""
    }
  }
  
  xmlQuery = paste(xmlQuery, linkedAttributeXML, linkedFilterXML,"</Dataset></Query>",sep="")
  
  if(verbose){
    cat(paste(xmlQuery,"\n", sep=""))
  }
  postRes = postForm(paste(martHost(mart),"?",sep=""),"query"=xmlQuery)
  
  if(postRes != ""){
    con = textConnection(postRes)
    result = read.table(con, sep="\t", header=bmHeader, quote = "\"", comment.char = "", as.is=TRUE, check.names = TRUE)
    close(con)
    if(all(is.na(result[,ncol(result)])))
      result = result[,-ncol(result),drop=FALSE]
  } else {
    warning("getLDS returns NULL.")
    result=NULL
  }
  return(result)
} 

######################
#getXML
######################

getXML <- function(host="http://www.biomart.org/biomart/martservice?", xmlquery){
      pf = postForm(host,"query"=xmlquery)
      con = textConnection(pf)
      result = read.table(con, sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
      close(con)
      return(result)
}

######################
#getBMlist
######################

getBMlist <- function(attributes, filters = "", values = "", mart, list.names = NULL, na.value = NA, verbose=FALSE, giveWarning=TRUE){
  if(giveWarning) writeLines("Performing your query using getBM is preferred as getBMlist perfoms a separate getBM query for each of the values one gives.  This is ok for a short list but will definitely fail when used with longer lists.  Ideally one does a batch query with getBM and then iterates over that result.")
  out <- vector("list", length(attributes))
  if(is.null(list.names))
    names(out) <- attributes
  else
    names(out) <- list.names
  for(j in seq(along = attributes)){
    tmp2 <- vector("list", length(values))
    names(tmp2) <- values
    for(k in seq(along = tmp2)){
      tst <- getBM(attributes = attributes[j], filters=filters, values = values[k], mart = mart, verbose = verbose)
      
      if(class(tst) == "data.frame"){
        tmp <- unlist(unique(tst[!is.na(tst)]), use.names = FALSE)
        if(length(tmp) > 0)
          tmp2[[k]] <- tmp
        else
          tmp2[[k]] <- na.value
      }else{
        tmp2[[k]] <- na.value
      }
      out[[j]] <- tmp2
    }
  }
  return(out)
}


###############################
#                             #
#Ensembl specific functions   #
###############################

getGene <- function( id, type, mart){
  martCheck(mart,"ensembl") 
  checkWrapperArgs(id, type, mart)
  symbolAttrib = switch(strsplit(martDataset(mart), "_", fixed = TRUE, useBytes = TRUE)[[1]][1],hsapiens = "hgnc_symbol",mmusculus = "mgi_symbol","external_gene_id")
  typeAttrib = switch(type,affy_hg_u133a_2 = "affy_hg_u133a_v2",type)
  attrib = c(typeAttrib,symbolAttrib,"description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id")
  table = getBM(attributes = attrib,filters = type, values = id, mart=mart)
  return(table)
}

getSequence <- function(chromosome, start, end, id, type, seqType, upstream, downstream, mart, verbose=FALSE){
  martCheck(mart,"ensembl")
  if(missing(seqType) || !seqType %in% c("cdna","peptide","3utr","5utr", "gene_exon", "transcript_exon","transcript_exon_intron","gene_exon_intron","coding","coding_transcript_flank","coding_gene_flank","transcript_flank","gene_flank")){
    stop("Please specify the type of sequence that needs to be retrieved when using biomaRt in web service mode.  Choose either gene_exon, transcript_exon,transcript_exon_intron, gene_exon_intron, cdna, coding,coding_transcript_flank,coding_gene_flank,transcript_flank,gene_flank,peptide, 3utr or 5utr")
  }
  if(missing(type))stop("Please specify the type argument.  If you use chromosomal coordinates to retrieve sequences, then the type argument will specify the type of gene indentifiers that you will retrieve with the sequences.  If you use a vector of identifiers to retrieve the sequences, the type argument specifies the type of identifiers you are using.")
  if(missing(id) && missing(chromosome) && !missing(type))stop("No vector of identifiers given. Please use the id argument to give a vector of identifiers for which you want to retrieve the sequences.")
  if(!missing(chromosome) && !missing(id))stop("The getSequence function retrieves sequences given a vector of identifiers specified with the id argument of a type specified by the type argument.  Or alternatively getSequence retrieves sequences given a chromosome, a start and a stop position on the chromosome.  As you specified both a vector of identifiers and chromsomal coordinates. Your query won't be processed.")
    
  if(!missing(chromosome)){
    if(!missing(start) && missing(end))stop("You specified a chromosomal start position but no end position.  Please also specify a chromosomal end position.")
    if(!missing(end) && missing(start))stop("You specified a chromosomal end position but no start position.  Please also specify a chromosomal start position.")
    if(!missing(start)){ start = as.integer(start)
                           end = as.integer(end)
                       }
    if(missing(upstream) && missing(downstream)){
      sequence = getBM(c(seqType,type), filters = c("chromosome_name","start","end"), values = list(chromosome, start, end), mart = mart, checkFilters = FALSE, verbose=verbose)
        }
    else{
      if(!missing(upstream) && missing(downstream)){
        sequence = getBM(c(seqType,type), filters = c("chromosome_name","start","end","upstream_flank"), values = list(chromosome, start, end, upstream), mart = mart, checkFilters = FALSE, verbose=verbose)
      }
      if(!missing(downstream) && missing(upstream)){
        sequence = getBM(c(seqType,type), filters = c("chromosome_name","start","end","downstream_flank"), values = list(chromosome, start, end, downstream), mart = mart, checkFilters = FALSE, verbose = verbose)
      }
      if(!missing(downstream) && !missing(upstream)){
        stop("Currently getSequence only allows the user to specify either an upstream of a downstream argument but not both.")
      }
    }
  }

  if(!missing(id)){
    if(missing(type)) stop("Type argument is missing.  This will be used to retrieve an identifier along with the sequence so one knows which gene it is from.  Use the listFilters function to select a valid type argument.")
    if(!type %in% listFilters(mart, what="name")) stop("Invalid type argument.  Use the listFilters function to select a valid type argument.")
  
    valuesString = paste(id,"",collapse=",",sep="")
    if(missing(upstream) && missing(downstream)){
      sequence = getBM(c(seqType,type), filters = type, values = id, mart = mart, verbose=verbose)
    }
    else{
      if(!missing(upstream) && missing(downstream)){
        sequence = getBM(c(seqType,type), filters = c(type, "upstream_flank"), values = list(id, upstream), mart = mart, checkFilters = FALSE, verbose=verbose)
      }
      if(!missing(downstream) && missing(upstream)){
        sequence = getBM(c(seqType,type), filters = c(type, "downstream_flank"), values = list(id, downstream), mart = mart, checkFilters = FALSE, verbose=verbose)
      }
      if(!missing(downstream) && !missing(upstream)){
        stop("Currently getSequence only allows the user to specify either an upstream of a downstream argument but not both.")
      }
    }
  }
  return(sequence)
}

####################
#export FASTA      #
####################

exportFASTA <- function( sequences, file ){
  if( missing( sequences ) || class( sequences ) != "data.frame"){
    stop("No data.frame given to write FASTA.  The data.frame should be the output of the getSequence function.");
  }
  if( missing(file)){
    stop("Please provide filename to write to");
  }
  if(length(sequences[1,]) == 2){
   for(i in seq(along = sequences[,2])){
     cat(paste(">",sequences[i,2],"\n",sep=""),file = file, append=TRUE);
     cat(as.character(sequences[i,1]),file = file, append = TRUE);
     cat("\n\n", file = file, append = TRUE);
   }
  }
  else{
   for(i in seq(along = sequences[,2])){
     cat(paste(">chromosome_",sequences[i,1],"_start_",sequences[i,2],"_end_",sequences[i,3],"\n",sep=""),file = file, append=TRUE);
     cat(as.character(sequences[i,4]),file = file, append = TRUE);
     cat("\n\n", file = file, append = TRUE);
   }
 }  
}

##########################
#Old function stubs
##########################

###################
#Nature Protocol
###################

NP2009code <- function(){
  edit(file=system.file('scripts', 'Integration-NP.R', package = 'biomaRt'))
}
