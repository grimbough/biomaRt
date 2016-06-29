##########################################
#getBM: generic BioMart query function   # 
##########################################

getBM <- function(attributes, filters = "", values = "", mart, curl = NULL, checkFilters = TRUE, verbose=FALSE, uniqueRows=TRUE, bmHeader=FALSE, quote="\""){
    
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
    
    xmlQuery = paste("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = '",martVSchema(mart),"' uniqueRows = '",as.numeric(uniqueRows),"' count = '0' datasetConfigVersion = '0.6' header='",as.numeric(bmHeader),"' requestid= 'biomaRt'> <Dataset name = '",martDataset(mart),"'>",sep="")
    
    #checking the Attributes
    invalid = !(attributes %in% listAttributes(mart, what="name"))
    if(any(invalid))
        stop(paste("Invalid attribute(s):", paste(attributes[invalid], collapse=", "),
                   "\nPlease use the function 'listAttributes' to get valid attribute names"))
    
    #check if attributes come from multiple attribute pages currently disabled until ID issue resovled at Ensembl
    # if(FALSE){
    #     att = listAttributes(mart, what=c("name","page"))
    #     att = att[which(att[,1] %in% attributes),]
    #     attOK = FALSE
    #     pages = unique(att[,2])
    #     if(length(pages) <= 1){
    #         attOK = TRUE
    #     }
    #     else{
    #         for(page in pages){
    #             if(length(attributes) == length(which(attributes %in% att[which(att[,2] == page),1]))) attOK = TRUE
    #         }
    #     }
    #     if(!attOK){
    #         stop(paste("Querying attributes from multiple attribute pages is not allowed.  To see the attribute pages attributes belong to, use the function attributePages."))
    #     }
    # }
    #attribute are ok lets add them to the query
    attributeXML =  paste("<Attribute name = '", attributes, "'/>", collapse="", sep="")
    
    #checking the filters
    if(filters[1] != "" && checkFilters){
        invalid = !(filters %in% listFilters(mart, what="name"))
        if(any(invalid))
            stop(paste("Invalid filters(s):", paste(filters[invalid], collapse=", "),
                       "\nPlease use the function 'listFilters' to get valid filter names"))
    }
    filterXML <- .generateFilterXML(filters, values, mart)
    
    
    xmlQuery = paste(xmlQuery, attributeXML, filterXML,"</Dataset></Query>",sep="")
    
    if(verbose){
        cat(paste(xmlQuery,"\n", sep=""))
    }      
    
    postRes = tryCatch(postForm(paste(martHost(mart),"?",sep=""),"query" = xmlQuery), error = function(e) {
        stop("Request to BioMart web service failed. Verify if you are still connected to the internet.  Alternatively the BioMart web service is temporarily down.")
        })
    
    if(verbose){
        writeLines("#################\nResults from server:")
        print(postRes)
    }
    if(!(is.character(postRes) && (length(postRes)==1L)))
        stop("The query to the BioMart webservice returned an invalid result: biomaRt expected a character string of length 1. Please report this to the mailing list.")
    
    if(gsub("\n", "", postRes, fixed = TRUE, useBytes = TRUE) == "") { # meaning an empty result
        
        result = as.data.frame(matrix("", ncol=length(attributes), nrow=0), stringsAsFactors=FALSE)
        
    } else {
        
        if(length(grep("^Query ERROR", postRes))>0L)
            stop(postRes)
        
        ## convert the serialized table into a dataframe
        con = textConnection(postRes)
        result = read.table(con, sep="\t", header=bmHeader, quote = quote, comment.char = "", check.names = FALSE, stringsAsFactors=FALSE)
        if(verbose){
            writeLines("#################\nParsed results:")
            print(result)
        }
        close(con)
        
        if(!(is(result, "data.frame") && (ncol(result)==length(attributes)))) {
            print(head(result))
            stop("The query to the BioMart webservice returned an invalid result: the number of columns in the result table does not equal the number of attributes in the query. Please report this to the mailing list.")
        }
    }
    if(!bmHeader){  #assumes order of results same as order of attibutes in input 
        colnames(result) = attributes
    }
    else{
        result <- .setResultColNames(result = result, mart = mart)
    }
    return(result)
}

.setResultColNames <- function(result, mart) {
    
    att = listAttributes(mart)
    resultNames = colnames(result)
    
    matches <- match(resultNames, att[,2], NA)
    if(any(is.na(matches))) {
        warning("Problems assigning column names. Currently using the biomart description field.  You may wish to set these manually.")
    } else {
        colnames(result) = att[matches, 1]
    }
    return(result)
}

.generateFilterXML <- function(filters = "", values, mart) {
    
    ## return emptry string if no filter specified
    if(filters[1]== "") {
        return("")
    }
    ## if we have multiple filters, the values must be specified as a list.
    if(length(filters) > 1 && class(values) != "list") {
        stop("If using multiple filters, the 'value' has to be a list.\nFor example, a valid list for 'value' could be: list(affyid=c('1939_at','1000_at'), chromosome= '16')\nHere we select on Affymetrix identifier and chromosome, only results that pass both filters will be returned");
    } 
    
    if(!is.list(values)){
        values <- list(values)
    }
    names(values) <- filters

    individualFilters <- sapply(filters, 
           function(filter, values, mart) {
                ## if the filter exists and is boolean we do this
                if(filter %in% listFilters(mart, what = "name") && grepl('boolean', filterType(filter = filter, mart = mart)) ) {
                    if(!is.logical(values[[filter]])) 
                            stop("biomaRt error: ", filter, " is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved.")
                        val <- ifelse(values[[filter]], yes = 0, no = 1)
                        val <- paste0("' excluded = \"", val, "\" ")
                } else { ## otherwise the filter isn't boolean, or doesn't exist, and we treat them the same
                    ## convert floats to integers
                    if(is.numeric(values[[filter]])) 
                        values[[filter]] <- as.integer(values[[filter]])
                    val <- paste0(values[[filter]], collapse = ",")
                    val <- paste0("' value = '", val, "' ")
                }
                filterXML <- paste0("<Filter name = '", filter, val, "/>")
                return(filterXML)
            }, values, mart)
    
    filterXML <- paste0(individualFilters, collapse = "")
    return(filterXML)
}

