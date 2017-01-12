
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