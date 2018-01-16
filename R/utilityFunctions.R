
## sometimes results can be returned by getMB() in a different order to we 
## asked for them, which messes up the column names.  Here we try to match
## results to known attribute names and rename accordingly.
.setResultColNames <- function(result, mart, attributes, bmHeader = FALSE) {
    
    ## get all avaialble sttributes and 
    ## filter only for the ones we've actually asked for
    att <- listAttributes(mart, what = c("name", "description"))
    att <- att[which(att[,'name'] %in% attributes), ]
    if(length(which(duplicated(att[,'description']))) > 
       length(which(duplicated(att)))) {
        warning("Cannot unambiguously match attribute names
                Ignoring bmHeader argument and using biomart 
                description field")
        return(result)
    }
    
    resultNames = colnames(result)
    ## match the returned column names with the attribute names
    matches <- match(resultNames, att[,2], NA)
    if(any(is.na(matches))) {
        warning("Problems assigning column names.",
                "Currently using the biomart description field.", 
                "You may wish to set these manually.")
        return(result)
    }
    ## if we want to use the attribute names we specified, do this, 
    ## otherwise we use the header returned with the query
    if(!bmHeader) {
        colnames(result) = att[matches, 1]
    }
    ## now put things in the order we actually asked for the attributes in
    result <- result[, match(att[matches,1], attributes), drop=FALSE]
    
    return(result)
}

## BioMart doesn't work well if the list of values provided to a filter is 
## longer than 500 values.  It returns only a subset of the requested data
## and does so silently!  This function is designed to take a list of provided
## filters, and split any longer than 'maxChunkSize'.  It operates recursively
## incase there are multiple filters that need splitting, and should ensure
## all possible groupings of filters are retained.
.splitValues <- function(valuesList, maxChunkSize = 500) {
    
    vLength <- vapply(valuesList[[1]], FUN = length, FUN.VALUE = integer(1))
    
    if(all(vLength <= maxChunkSize)) {
        return(valuesList)
    } else {
        ## pick the next filter to split
        vIdx <- min(which(vLength > maxChunkSize))
        
        nchunks <- (vLength[vIdx] %/% maxChunkSize) + 1
        splitIdx <- rep(1:nchunks, each = ceiling(vLength[vIdx] / nchunks))[ 1:vLength[vIdx] ]
        
        ## a new list we will populate with the chunks
        tmpList <- list()
        for(i in 1:nchunks) {
            for( j in 1:length(valuesList) ) {
                listIdx <- ((i - 1) * length(valuesList)) + j
                tmpList[[ listIdx ]] <- valuesList[[j]]
                tmpList[[ listIdx ]][[ vIdx ]] <- tmpList[[ listIdx ]][[ vIdx ]][which(splitIdx == i)]
            }
        }
        ## recursively call the function to process next filter
        valuesList <- .splitValues(tmpList)
    }
    return(valuesList)
}

## Creating the filter XML for a single chunk of values.  Returns a character
## vector containing the XML lines for all specified filters & their 
## attributes spliced together into a single string.
.createFilterXMLchunk <- function(filterChunk, mart) {
    
    individualFilters <- vapply(names(filterChunk), 
        FUN = function(filter, values, mart) {
            
            ## if the filter exists and is boolean we do this
            if(filter %in% listFilters(mart, what = "name") && 
               grepl('boolean', filterType(filter = filter, mart = mart)) ) {
                if(!is.logical(values[[filter]])) 
                    stop("biomaRt error:\n", 
                         filter, " is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved.")
                val <- ifelse(values[[filter]], yes = 0, no = 1)
                val <- paste0("' excluded = \"", val, "\" ")
                
            } else { 
                ## otherwise the filter isn't boolean, or doesn't exist
                
                if(is.numeric(values[[filter]])) 
                    values[[filter]] <- as.integer(values[[filter]])
                val <- paste0(values[[filter]], collapse = ",")
                val <- paste0("' value = '", val, "' ")
            }
            filterXML <- paste0("<Filter name = '", filter, val, "/>")
            return(filterXML)
        }, FUN.VALUE = character(1), 
        filterChunk, mart,
        USE.NAMES = FALSE)
    
    filterXML <- paste0(individualFilters, collapse = "")
    return(filterXML)
}

.generateFilterXML <- function(filters = "", values, mart) {
    
    ## return empty string if no filter specified & this isn't ensembl
    ## specifying no filter is generally bad, as it will get 'everything'
    ## and we might encounter the time out problem
    if(filters[1] == "") {
        return("")
    }
    ## if we have multiple filters, the values must be specified as a list.
    if(length(filters) > 1 && class(values) != "list") {
        stop("If using multiple filters, the 'value' has to be a list.\nFor example, a valid list for 'value' could be: list(affyid=c('1939_at','1000_at'), chromosome= '16')\nHere we select on Affymetrix identifier and chromosome, only results that pass both filters will be returned");
    } 
    ## it's easy to not realise you're passing a data frame here, so check
    if(is.data.frame(values) && ncol(values == 1)) {
        values <- values[,1]
    }

    
    if(!is.list(values)){
        values <- list(values)
    }
    names(values) <- filters
    
    values <- .splitValues(list(values))
    
    filterXML_list <- lapply(values, .createFilterXMLchunk, mart)
}

#' it seems like pretty common practice for users to copy and paste the host
#' name from a browser if they're not accessing Ensembl.  Typically this will
#' include the "http://" and maybe a trailing "/" and this messes up our
#' paste the complete URL strategy and produces something invalid.  
#' This function tidies that up to catch common variants.
.cleanHostURL <- function(host) {
    
    ## strip trailing slash
    host <- gsub(pattern = "/$", replacement = "", x = host)
    
    ## just supplying 'ensembl.org' is no longer handled correctly
    ## stick 'www' infront if we see this
    if( grepl(pattern = "^ensembl\\.org$", x = host) ) {
        host = "www.ensembl.org"
    }
    
    ## only prepend http if needed 
    if(!grepl(pattern = "^http://|^https://", x = host)) {
        host <- paste0("http://", host)
    }
    
    return(host)
}

#' ensembl redirection doesn't seem to be working properly as of 12-12-2017
#' This is a wrapper function to catch POSTS that are redirected and fail
#' The new host is captured from the header and used in a re-submission
.submitQuery <- function(host, query) {
    res <- httr::POST(url = host,
                      body = list('query' = query),
                      config = list(timeout(30)))

    if(res$all_headers[[1]]$status == 302) {
        host <- stringr::str_match(string = res$all_headers[[1]]$headers$location,
                               pattern = "//([a-zA-Z./]+)\\??;?redirectsrc")[,2]
        res <- httr::POST(url = host,
                          body = list('query' = query),
                          config = list(timeout(30)))
    }
    ## content() prints a message about encoding not being supplied 
    ## for ensembl.org - no default, so we suppress it
    return( suppressMessages(content(res)) )
}


