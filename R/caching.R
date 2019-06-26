###########################################################
## Functions for caching
###########################################################

.createHash <- function(mart, attributes, filters, values) {
    
    ## if we are using the current version Ensembl URL
    ## swap for the archive version so we can check when it is outdated
    host <- mart@host
    if(grepl("(www|uswest|useast|asia)\\.ensembl\\.org", host)) {
        archives <- biomaRt::listEnsemblArchives()
        host <- archives[which(archives$current_release == "*"), "url"]
    }
    
    attributes <- paste( sort(attributes), collapse = "" )
    idx <- order(filters)
    filters <- paste( filters[idx], collapse = "" )
    if(is.list(values)) {
        values <- values[idx]
        values <- unlist(lapply(values, sort))
    } else {
        values <- sort(values)
    }
    values <- paste( values, collapse = "" )    
    
    combined <- paste(c(host, mart@biomart, mart@dataset, attributes, filters, values), 
                      collapse = "_")
    paste0("biomaRt_", 
           as(openssl::md5(combined), "character"))
}


.checkCache <- function(bfc, hash) {
    res <- bfcquery(bfc, query = hash, field = "rname")
    as.logical(nrow(res))
}

clearBiomartCache <- function() {
    cache <- rappdirs::user_cache_dir(appname="biomaRt")
    bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
    removebfc(bfc, ask = FALSE)
}
