###########################################################
## Functions for caching
###########################################################

createHash <- function(mart, attributes, filters, values) {
    
    ## if we are using the current version Ensembl URL
    ## swap for the archive version so we can check when it is outdated
    host <- mart@host
    if(grepl("(www|uswest|useast|asia)\\.ensembl\\.org", host)) {
        archives <- biomaRt::listEnsemblArchives()
        host <- archives[which(archives$current_release == "*"), "url"]
    }
    
    tmp <- lapply( list(attributes, filters, values), 
                   FUN = function(x) {
                       if(is.list(x))
                           x <- unlist(lapply(x, sort))
                       else
                           x <- sort(x)
                       paste( x, collapse = "")
                   })
    
    combined <- paste(c(host, unlist(tmp)), 
                      collapse = "")
    as(openssl::md5(combined), "character")
}


checkCache <- function(bfc, hash) {
    res <- bfcquery(bfc, query = hash, field = "rname")
    as.logical(nrow(res))
}
