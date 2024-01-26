
##########################
#biomaRt source code     #
##########################
#                        #
#Licence: Artistic       #
#Author: Steffen Durinck #
##########################


##############################################################
#martCheck                                                   # 
#                                                            #
#This function checks if there is a valid Mart object,       # 
#if a dataset is selected and                                #
#if the correct BioMart database has been selected (optional)# 
##############################################################

martCheck = function(mart, biomart = NULL){
    if( missing(mart) || !inherits(mart,"Mart") ) {
        stop("You must provide a valid Mart object. To create a Mart object use the function: useMart.  Check ?useMart for more information.")
    }
    if(!is.null(biomart)){
        martcheck = martBM(mart)
        bmok = FALSE
        for(k in seq_along(biomart)) {
            if(martcheck[1] == biomart[k]) {	
                bmok = TRUE
            }
        }		    
        if(!bmok){
            stop(paste("This function only works when used with the ",biomart," BioMart.",sep="")) 
        }      
    }
    if(martDataset(mart)==""){
        stop("No dataset selected, please select a dataset first.  You can see the available datasets by using the listDatasets function see ?listDatasets for more information.  Then you should create the Mart object by using the useMart function.  See ?useMart for more information");
    }
}


bmRequest <- function(request, httr_config, verbose = FALSE){
    if(verbose) 
        message("Attempting web service request:\n", request)
    
    result <- httr::GET(request, config = httr_config, 
                        content_type("text/plain"),
                        timeout(getOption("timeout", default = 60)))
    stop_for_status(result)

    result2 <- content(result, encoding = "UTF-8", as = "text")
    if(is.na(result2)) {
        result2 <- content(result, encoding = "Latin1", as = "text")
    }
    return(result2)
}

#######################################################
#listMarts:                                           #
#list all available BioMart databases by default      #
#listMarts will check the central service to see which#
#BioMart databases are present                        #
#######################################################

listMarts <- function( mart = NULL, host="https://www.ensembl.org", path="/biomart/martservice", 
                       port, includeHosts = FALSE, archive = FALSE, httr_config, verbose = FALSE){
    
    if(missing(port)) {
        port <- ifelse(grepl("https", host), yes = 443, no = 80)
    }
    
    if(grepl(pattern = "^https://.*ensembl.org", x = host) && missing(httr_config)) {
      httr_config <- .getEnsemblSSL()
    }
  
    if(missing(httr_config)) {
        httr_config <- httr::config()
    }
    
    .listMarts(mart = mart, host = host, path = path, port = port, includeHosts = includeHosts,
                archive = archive, verbose = verbose, httr_config = httr_config, ensemblRedirect = TRUE)
    
}

.listMarts <- function( mart = NULL, host="www.ensembl.org", path="/biomart/martservice", 
                       port=80, includeHosts = FALSE, archive = FALSE, verbose = FALSE, 
                       httr_config, ensemblRedirect = NULL, warn = TRUE){

    request = NULL
    if(is.null(mart)){
        host <- .cleanHostURL(host, warn = warn)
        if(archive) {
            stop("The archive = TRUE argument is now defunct.\n", 
                 "Use listEnsemblArchives() to find the URL to directly query an Ensembl archive.")
        } else {
            request <- paste0(host, ":", port, path, "?type=registry&requestid=biomaRt")
        }
        if(is(httr_config, 'list')) {
            httr_config <- do.call(c, httr_config)
        }
    } else if(is(mart, 'Mart')) {
            request = paste0(martHost(mart), "?type=registry&requestid=biomaRt") 
            httr_config <- martHTTRConfig(mart)
    } else{
            stop(mart, " object needs to be of class Mart created with the useMart function.\n",
            "If you don't have a Mart object yet, use listMarts() without arguments or only specify the host argument")
    } 	
    
    if(!ensemblRedirect && grepl(x = request, pattern = "ensembl.org")) {
        request <- paste0(request, "&redirect=no")
    }
    
    registry = bmRequest(request = request, httr_config = httr_config, verbose = verbose)
    
    ## check this looks like the MartRegistry XML, otherwise throw an error
    if(!grepl(x = registry, pattern = "^\n*<MartRegistry>")) {
        
        if(grepl(x = registry, pattern = "status.ensembl.org")) {
            stop("Your query has been redirected to https://status.ensembl.org ",
                 "indicating this Ensembl service is currently unavailable.",
                 "\nLook at ?useEnsembl for details on how to try a mirror site.",
                 call. = FALSE)
        } else {
            stop('Unexpected format to the list of available marts.\n',
                 'Please check the following URL manually, ',
                 'and try ?listMarts for advice.\n',
                 request, 
                 call. = FALSE)
        }
    }
    registry = xmlTreeParse(registry, asText=TRUE)
    registry = registry$doc$children[[1]]
    
    marts = list(biomart = NULL, version = NULL, host = NULL, path = NULL, database = NULL)
    index = 1
    
    # if(host != "www.biomart.org" || archive){
        for(i in seq(length.out=xmlSize(registry))){
            if(xmlName(registry[[i]])=="MartURLLocation"){  
                if(xmlGetAttr(registry[[i]],"visible") == 1){
                    if(!is.null(xmlGetAttr(registry[[i]],"name"))) marts$biomart[index] = as.character(xmlGetAttr(registry[[i]],"name"))
                    if(!is.null(xmlGetAttr(registry[[i]],"database"))) marts$database[index] = as.character(xmlGetAttr(registry[[i]],"database"))
                    if(!is.null(xmlGetAttr(registry[[i]],"displayName"))) marts$version[index] = as.character(xmlGetAttr(registry[[i]],"displayName"))
                    if(!is.null(xmlGetAttr(registry[[i]],"host"))) marts$host[index] = as.character(xmlGetAttr(registry[[i]],"host"))
                    if(!is.null(xmlGetAttr(registry[[i]],"path"))) marts$path[index] = as.character(xmlGetAttr(registry[[i]],"path"))
                    if(!is.null(xmlGetAttr(registry[[i]],"port"))) marts$port[index] = as.character(xmlGetAttr(registry[[i]],"port"))
                    if(!is.null(xmlGetAttr(registry[[i]],"serverVirtualSchema"))){
                        marts$vschema[index] =  as.character(xmlGetAttr(registry[[i]],"serverVirtualSchema"))
                    }
                    index=index+1
                }
            }
        }
    if(includeHosts){
        return(marts)
    }
    else{
        ret = data.frame(biomart = as.character(marts$biomart),
                         version = as.character(marts$version), 
                         stringsAsFactors=FALSE)
        return(ret)
    } 
}

#################################
# #                           # #
# # Generic BioMart functions # #
# #                           # #
#################################

useMart <- function(biomart, dataset, host = "https://www.ensembl.org", path = "/biomart/martservice", port, 
                     archive = FALSE, version, verbose = FALSE) {
    
    if(missing(port)) {
        port <- ifelse(grepl("https", host)[1], yes = 443, no = 80)
    }
    
    mart <- .useMart(biomart, dataset, host = host, path = path, port = port, 
                     archive = archive, version = version, verbose = verbose, 
                     httr_config = list(httr::config()), ensemblRedirect = TRUE)
}

.useMart <- function(biomart, dataset, host = "https://www.ensembl.org", path = "/biomart/martservice", port = 443, 
                    archive = FALSE, ensemblRedirect = NULL, version, httr_config, verbose = FALSE){
    
    if(missing(biomart) && missing(version)) 
        stop("No biomart databases specified. Specify a biomart database to use using the biomart or version argument")
    if(!missing(biomart)){ 
        if(!(is.character(biomart)))
            stop("biomart argument is not a string. ",
                 "The biomart argument should be a single character string")
    }
    
    if(biomart == "ensembl" & grepl(x = host, pattern = "ensembl.org")) {
        biomart = "ENSEMBL_MART_ENSEMBL"
    }
    
    reqHost = host
    host <- .cleanHostURL(host)
    
    marts <- .listMarts(host=host, path=path, port=port, includeHosts = TRUE,
                       httr_config = httr_config, archive = archive,
                       ensemblRedirect = ensemblRedirect, warn = FALSE)
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
    
    if(is.na(marts$path[mindex]) || is.na(marts$vschema[mindex]) || 
       is.na(marts$host[mindex]) || is.na(marts$port[mindex])) 
        stop("The selected biomart databases is not available due to error in the BioMart central registry, please report so the BioMart registry file can be fixed.")
    
    if(marts$path[mindex]=="") marts$path[mindex]="/biomart/martservice" #temporary to catch bugs in registry

    if(!missing(version)) biomart = marts$biomart[mindex]
    biomart = sub(" ","%20",biomart, fixed = TRUE, useBytes = TRUE)
    
    ## adding option to force use of specified host with ensembl
    redirect <- ifelse(!ensemblRedirect && grepl(x = host, pattern = "ensembl.org"), 
                       "?redirect=no",
                       "")
    
    if(missing(httr_config)) {
        httr_config <- list()
    }
    
    mart <- Mart( 
                biomart = biomart,
                vschema = marts$vschema[mindex], 
                host = paste0(host, ":", 
                              port,
                              marts$path[mindex],
                              redirect),
                httr_config = httr_config
            )
    
    if(length(grep("archive",martHost(mart)) > 0)){
        
        ## hack to work around redirection of most recent mirror URL
        archives <- .listEnsemblArchives(https = TRUE, httr_config = httr_config)
        current_release <- archives[archives$current_release == "*", 'url']
        if(grepl(martHost(mart), pattern = current_release)) {
            martHost(mart) <- stringr::str_replace(martHost(mart), pattern = current_release, "https://www.ensembl.org")
            martHost(mart) <- stringr::str_replace(martHost(mart), pattern = ":80/", ":443/")
        }
    }
    
    BioMartVersion=bmVersion(mart, verbose=verbose)
    
    if(verbose){
        writeLines(paste("BioMartServer running BioMart version:",BioMartVersion,sep=" "))
        writeLines(paste("Mart virtual schema:",martVSchema(mart),sep=" "))
        if(length(grep(reqHost,martHost(mart))) == 0){
            writeLines(paste("Requested host was redirected from ", reqHost, " to " ,martHost(mart),sep=""))
        } 
        writeLines(paste("Mart host:",martHost(mart),sep=" "))
    }
    if(!missing(dataset)){
        mart = useDataset(mart = mart, dataset=dataset, verbose = verbose)
    }
    return(mart)
}

listDatasets <- function(mart, verbose = FALSE) {
    .listDatasets(mart = mart, verbose = verbose, sort = TRUE)
}

.listDatasets <- function(mart, verbose = FALSE, sort = FALSE) {
    if(missing(mart) || !is(mart, 'Mart'))
        stop("No Mart object given or object not of class 'Mart'")
    
    ## we choose a separator based on whether 'redirect=no' is present
    ## should always be '?' now
    sep <- ifelse(grepl(x = martHost(mart), pattern = ".+\\?.+"), "&", "?")
    
    request = paste0(martHost(mart), sep, "type=datasets&requestid=biomaRt&mart=", martBM(mart))
    httr_config <- martHTTRConfig(mart)
    
    bmResult = bmRequest(request = request, httr_config = httr_config, verbose = verbose)
    con = textConnection(bmResult)
    txt = scan(con, sep="\t", blank.lines.skip=TRUE, what="character", quiet=TRUE, quote = "\"")
    close(con)
    
    ## select visible ("1") table sets
    i = intersect(which(txt=="TableSet"), which(txt=="1")-3L)
    
    res = data.frame(dataset     = I(txt[i+1L]),
                     description = I(txt[i+2L]),
                     version     = I(txt[i+4L]))
    
    ## sort alphabetically
    if(sort)
        res <- res[ order(res$dataset), ]
    rownames(res) <- NULL
    
    return(res)
}

## Check version of BioMart service
bmVersion <- function(mart, verbose=FALSE){
    
    ## save some time and a HTTP request if this is Ensembl
    if(grepl(pattern = "ensembl.org", x = martHost(mart), fixed = TRUE)) {
        bmv <- "0.7"
    } else {
        ## we choose a separator based on whether 'redirect=no' is present
        sep <- ifelse(grepl(x = martHost(mart), pattern = ".+\\?.+"), "&", "?")
        
        request = paste0(martHost(mart), sep, "type=version", "&requestid=biomaRt&mart=", martBM(mart))
        httr_config <- martHTTRConfig(mart)
        
        BioMartVersion = bmRequest(request = request, httr_config = httr_config, verbose = verbose)
        bmv = ""
        if(BioMartVersion == "\n" || BioMartVersion == ""){
            bmv = NA
            if(verbose) warning(paste("BioMart version is not available from BioMart server:",request,sep="\n"))
        }
        else{
            con = textConnection(BioMartVersion)
            bmVersionParsed = read.table(con, sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
            close(con)
            if(verbose) print(bmVersionParsed)
            
            if(dim(bmVersionParsed)[2] >= 1){
                bmv=bmVersionParsed[1,1]
            }
        }
    }
    return(bmv)
}


.getAttrFilt <- function(mart, verbose, type) {
    
    ## we choose a separator based on whether 'redirect=no' is present
    sep <- ifelse(grepl(x = mart@host, pattern = ".+\\?.+"), "&", "?")
    
    request <- paste0(mart@host, sep, "type=", type,
                     "&dataset=", martDataset(mart),
                     "&requestid=biomaRt&mart=", martBM(mart),
                     "&virtualSchema=", martVSchema(mart))

    attrfilt <- bmRequest(request = request, httr_config = martHTTRConfig(mart), verbose = verbose)
    attrfiltParsed <- read.table(text = attrfilt, sep="\t", header=FALSE, 
                                quote = "", comment.char = "", as.is=TRUE)
    return(attrfiltParsed)

}

.getAttributes <- function(mart, verbose = FALSE) {
    
    attributes_table <- .getAttrFilt(mart = mart, verbose = verbose, type = "attributes")
    
    if(ncol(attributes_table) < 4)
        stop("biomaRt error: looks like we're connecting to incompatible version of BioMart.")
    
    colnames(attributes_table) <- c("name", "description",
                                    "fullDescription", "page")
    return(attributes_table)
}

.getFilters <- function(mart, verbose = FALSE) {
    
    filters_table <- .getAttrFilt(mart = mart, verbose = verbose, type = "filters")
    
    if(ncol(filters_table) < 7)
        stop("biomaRt error: looks like we're connecting to incompatible version of BioMart.")
    
    colnames(filters_table) <- c("name", "description", "options",
                                 "fullDescription", "filters",
                                 "type", "operation")
    return(filters_table)
}

## Utility function to check dataset specification
## Returns dataset name as a character assuming all checks
## have been passed.
checkDataset <- function(dataset, mart) {
    
    validDatasets <- .listDatasets(mart, sort = FALSE)
    ## subsetting data.frames can produce some weird classes
    ## which aren't character(), so we coerce it here
    dataset <- as.character(dataset)
    
    if(length(dataset) > 1) 
        stop("Please only specify a single dataset name")
    
    if(is.na(match(dataset, validDatasets$dataset)))
        stop(paste("The given dataset: ",dataset,", is not valid.  Correct dataset names can be obtained with the listDatasets() function."))
    
    return(dataset)
}

## Select a BioMart dataset             
useDataset <- function(dataset, mart, verbose = FALSE){
    if( missing(mart) || !inherits(mart,"Mart") )
        stop("No valid Mart object given, specify a Mart object with the attribute mart")
    
    if(missing(dataset)) {
        stop("No dataset given.  Please use the dataset argument to specify which dataset you want to use. Correct dataset names can be obtained with the listDatasets() function.")
    } else {
        dataset <- checkDataset(dataset = dataset, mart = mart)
    }
    martDataset(mart) <- dataset  
    
    if(verbose) message("Checking attributes ...", appendLF = FALSE)
    martAttributes(mart) <- .getAttributes(mart, verbose = verbose)
    if(verbose){
        message(" ok")
        message("Checking filters ...", appendLF = FALSE)
    }
    martFilters(mart) <- .getFilters(mart, verbose = verbose)
    if(verbose) message(" ok")
    return( mart )
}

## listAttributes
listAttributes <- function(mart, page, what = c("name","description","page")) {
    martCheck(mart)
    if(!missing(page) && !page %in% attributePages(mart)) 
        stop("The chosen page: ",page," is not valid, please use the correct page name using the attributePages function")
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
listFilters <- function(mart, what = c("name", "description")) {
    
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
    .Deprecated(new = "listFilterOptions",
                msg = c("filterOptions() has been deprecated and will be removed from biomaRt.",
                "\nPlease use listFilterOptions() instead."))
    listFilterOptions(mart = mart, filter = filter)
}

## filterType
filterType <- function(filter, mart){
    if(missing(filter)) 
        stop("No filter given. Please specify the filter for which you want to retrieve the filter type")
    if(!is.character(filter))
        stop("Filter argument should be of class character")
    martCheck(mart)
    type="unknown"
    sel = which(listFilters(mart, what="name") == filter)
    if(is.null(sel))
        stop(paste("Invalid filter",filter, sep=": "))
    type = listFilters(mart, what="type")[sel]
    return(type)
}

##########################################
#getBM: generic BioMart query function   # 
##########################################

getBM <- function(attributes, filters = "", values = "", mart, curl = NULL, 
                  checkFilters = TRUE, verbose=FALSE, uniqueRows=TRUE, bmHeader=FALSE, quote="\"",
                  useCache = TRUE){
    
    ## check the arguments are all valid
    martCheck(mart)
    if(missing( attributes ))
        stop("Argument 'attributes' must be specified.")
    
    if(is.list(filters) && !missing( values ))
        warning("Argument 'values' should not be used when argument 'filters' is a list and will be ignored.")
    if(is.list(filters) && is.null(names(filters)))
        stop("Argument 'filters' must be a named list when sent as a list.")
    if(!is.list(filters) && all(filters != "") && missing( values ))
        stop("Argument 'values' must be specified.")
    if(length(filters) > 0 && length(values) == 0)
        stop("Values argument contains no data.")
    if(is.list(filters)){
        values = filters
        filters = names(filters)
    }
    if(!is.logical(uniqueRows))
        stop("Argument 'uniqueRows' must be a logical value, so either TRUE or FALSE")
    
    ## determine if we should use the results cache
    if(useCache) {
        cache <- .biomartCacheLocation()
        bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
    }
    hash <- .createHash(mart, attributes, filters, values, uniqueRows, bmHeader)
    if( useCache && .checkValidCache(bfc, hash) ) {
        
        if(verbose) {
            message("Cache found")
        }
        result <- .readFromCache(bfc, hash)
        return(result)

    } else { 
    
        ## force the query to return the 'descriptive text' header names with the result
        ## we use these later to match and order attribute/column names    
        xmlQuery = paste0("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = '",
                          martVSchema(mart),
                          "' uniqueRows = '",
                          as.numeric(uniqueRows),
                          "' count='0' datasetConfigVersion='0.6' header='1'",
                          " formatter='TSV' requestid='biomaRt'> <Dataset name = '",
                          martDataset(mart),"'>")
        
        #checking the Attributes
        invalid = !(attributes %in% listAttributes(mart, what="name"))
        if(any(invalid))
            stop(paste("Invalid attribute(s):", paste(attributes[invalid], collapse=", "),
                       "\nPlease use the function 'listAttributes' to get valid attribute names"))
        
        #attribute are ok lets add them to the query
        attributeXML = paste0("<Attribute name = '", attributes, "'/>", collapse="")
        
        #checking the filters
        if(filters[1] != "" && checkFilters){
            invalid = !(filters %in% listFilters(mart, what="name"))
            if(any(invalid))
                stop(paste("Invalid filters(s):", paste(filters[invalid], collapse=", "),
                           "\nPlease use the function 'listFilters' to get valid filter names"))
        }
        
        ## filterXML is a list containing filters with reduced numbers of values
        ## to meet the 500 value limit in BioMart queries
        filterXmlList <- .generateFilterXML(filters, values, mart)
        
        resultList <- list()
        if(length(filterXmlList) > 1) {
            pb <- progress_bar$new(total = length(filterXmlList),
                                   width = options()$width - 10,
                                   format = "Batch submitting query [:bar] :percent eta: :eta")
            pb$tick(0)
            on.exit( pb$terminate() )
        }
    
        ## we submit a query for each chunk of the filter list
        for(i in seq_along(filterXmlList)) {
            
            if(i > 1) {
                pb$tick()
            }
            
            filterXML <- filterXmlList[[ i ]]
            fullXmlQuery <- paste(xmlQuery, attributeXML, filterXML,"</Dataset></Query>",sep="")
            
            if(verbose) {
                message(fullXmlQuery)
            }      
            
            ## we choose a separator based on whether '?redirect=no' is present
            sep <- ifelse(grepl(x = martHost(mart), pattern = ".+\\?.+"), "&", "?")
            
            ## create a unique name for this chunk & see if it has been run before
            chunk_hash <- digest::digest(paste(martHost(mart), fullXmlQuery), algo = "md5", serialize = FALSE)
            tf <- file.path(tempdir(), paste0("biomaRt_tmp_", chunk_hash, ".rds"))
            if(!file.exists(tf)) {
                postRes <- .submitQueryXML(host = paste0(martHost(mart), sep),
                                       query = fullXmlQuery,
                                       httr_config = martHTTRConfig(mart))
                result <- .processResults(postRes, mart = mart, hostURLsep = sep, fullXmlQuery = fullXmlQuery,
                                          quote = quote, numAttributes = length(attributes))
                saveRDS(result, file = tf)
            } else {
                result <- readRDS(tf)
            }
            resultList[[i]] <- .setResultColNames(result, mart = mart, 
                                                  attributes = attributes, bmHeader = bmHeader)
        }
        ## collate results
        result <- do.call('rbind', resultList)
    
        if(useCache) {
            .addToCache(bfc = bfc, result = result, hash = hash)
        }
        
        ## remove any temp chunk files
        file.remove( list.files(tempdir(), pattern = "^biomaRt.*rds$", full.names = TRUE) )
        return(result)
    }
}

###################################
#getLDS: Multiple dataset linking #
###################################

getLDS <- function(attributes, filters = "", values = "", mart, 
                   attributesL, filtersL = "", valuesL = "", martL, 
                   verbose = FALSE, uniqueRows = TRUE, bmHeader = TRUE) {
    
    martCheck(mart)
    martCheck(martL)
    
    if(martHost(mart) != martHost(martL)) {
        stop('Both datasets must be located on the same host.')
    }
    
    if(martBM(mart) != martBM(martL)) {
        stop('Both datasets must be located in the same Mart.\n',
             'You are trying to combine datasets in ', 
             martBM(mart), ' and ', martBM(martL))
    }
    
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
    
    xmlQuery = sprintf("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query virtualSchemaName = '%s' uniqueRows = '%s' count = '0' datasetConfigVersion = '0.6' header='%s' formatter = 'TSV' requestid= 'biomaRt'> <Dataset name = '%s'>",
                       martVSchema(mart) ,as.numeric(uniqueRows), as.numeric(bmHeader), martDataset(mart))
    
    attributeXML = paste("<Attribute name = '", attributes, "'/>", collapse="", sep="")

    ## ignore the chunk size here
    filterXML <- .generateFilterXML(filters = filters, values = values, 
    								mart = mart, maxChunkSize = Inf)
    
    xmlQuery = paste0(xmlQuery, attributeXML, filterXML,"</Dataset>")

    xmlQuery = paste0(xmlQuery, "<Dataset name = '",martDataset(martL),"' >")
    linkedAttributeXML =  paste("<Attribute name = '", attributesL, "'/>", collapse="", sep="")  
    linkedFilterXML <- .generateFilterXML(filters = filtersL, values = valuesL, 
    									  mart = mart, maxChunkSize = Inf)
        
    xmlQuery = paste0(xmlQuery, linkedAttributeXML, linkedFilterXML,"</Dataset></Query>")
    
    if(verbose){
        message(xmlQuery)
    }

    ## we choose a separator based on whether '?redirect=no' is present
    sep <- ifelse(grepl(x = martHost(mart), pattern = ".+\\?.+"), "&", "?")
    ## POST query
    postRes <- .submitQueryXML(host = paste0(martHost(mart), sep),
                               query = xmlQuery,
                               httr_config = martHTTRConfig(mart))
    
    if(length(grep("^Query ERROR", postRes))>0L)
        stop(postRes)  

    if(postRes != ""){
        con = textConnection(postRes)
        result = read.table(con, sep="\t", header=bmHeader, quote = "\"", comment.char = "", as.is=TRUE, check.names = TRUE)
        close(con)
        
        if(nrow(result) > 0 && all(is.na(result[,ncol(result)])))
            result = result[,-ncol(result),drop=FALSE]

        res_attributes <- c(attributes,attributesL)
        if(!(is(result, "data.frame") && (ncol(result)==length(res_attributes)))) {
            print(head(result))
            stop("The query to the BioMart webservice returned an invalid result: ", 
            "the number of columns in the result table does not equal the number of attributes in the query. \n",
            "Please report this on the support site at http://support.bioconductor.org")
        } 
        if(!bmHeader){  #assumes order of results same as order of attibutes in input
            colnames(result) = res_attributes

        }
    } else {
        warning("getLDS returns NULL.")
        result=NULL
    }
    return(result)
} 

######################
#getBMlist
######################

getBMlist <- function(attributes, filters = "", values = "", mart, list.names = NULL, 
                      na.value = NA, verbose=FALSE, giveWarning=TRUE){
    .Defunct(new = "getBM",
             msg = c("getBMlist() has been removed from biomaRt",
                     "\nPlease use getBM() instead")
    )
}



####################
#export FASTA      #
####################

exportFASTA <- function( sequences, file ) {
    if( missing( sequences ) || !is.data.frame( sequences )) {
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

###################
#Nature Protocol
###################

NP2009code <- function(){
    edit(file=system.file('scripts', 'Integration-NP.R', package = 'biomaRt'))
}
