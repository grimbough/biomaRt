###########################################################
## Functions for caching
###########################################################

.createHash <- function(mart, attributes, filters, values, uniqueRows = TRUE, bmHeader = FALSE) {
    
    ## if we are using the current version Ensembl URL
    ## swap for the archive version so we can check when it is outdated
    host <- martHost(mart)
    if(grepl("(www|useast|asia)\\.ensembl\\.org", host)) {
        archives <- .listEnsemblArchives(httr_config = martHTTRConfig(mart))
        host <- archives[which(archives$current_release == "*"), "url"]
    }
    
    attributes <- paste( sort(attributes), collapse = "" )
    ## need to keep the filters and values in the same order
    ## so create a single index for reordering both
    idx <- order(filters)
    filters <- paste( filters[idx], collapse = "" )
    if(is.list(values)) {
        values <- values[idx]
        values <- unlist(lapply(values, sort))
    } else {
        values <- sort(values)
    }
    values <- paste( values, collapse = "" )    
    
    combined <- paste(c(host, mart@biomart, mart@dataset, attributes, filters, values, uniqueRows, bmHeader), 
                      collapse = "_")
    paste0("biomaRt_", digest::digest(combined, algo = "md5", serialize = FALSE))
}

#' @param bfc Object of class BiocFileCache, created by a call to 
#' BiocFileCache::BiocFileCache()
#' @param hash unique hash representing a query.
.addToCache <- function(bfc, result, hash) {
  
  if(!dir.exists(.biomartCacheLocation()))
    dir.create(.biomartCacheLocation())
  
  ## write our file to the biomart cache location directly
  tf <- tempfile(tmpdir = .biomartCacheLocation())
  saveRDS(result, file = tf)
  
  ## check once more that there isn't an entry with this hash
  ## if its free add our new file
  ## if there's a clash don't add anything and tidy up
  if(!.checkInCache(bfc, hash = hash)) {
    bfcadd(bfc, rname = hash, fpath = tf, action = "asis")
    res <- TRUE
  } else {
    file.remove(tf)
    res <- FALSE
  }
  return(invisible(res))
}

#' @param bfc Object of class BiocFileCache, created by a call to 
#' BiocFileCache::BiocFileCache()
#' @param hash unique hash representing a query.
.readFromCache <- function(bfc, hash) {

    cache_hits <- bfcquery(bfc, hash, field = "rname")
    if(nrow(cache_hits) > 1) {
        stop("Multiple cache results found.",
             "\nPlease clear your cache by running biomartCacheClear()")
    } else {
        rid <- cache_hits$rid
        result <- readRDS( bfc[[ rid ]] )
        return(result)
    }
}

#' @param bfc Object of class BiocFileCache, created by a call to 
#' BiocFileCache::BiocFileCache()
#' @param hash unique hash representing a query.
#' 
#' This function returns TRUE if a record with the requested hash already 
#' exists in the file cache, otherwise returns FALSE.
#' @keywords Internal
.checkInCache <- function(bfc, hash, verbose = FALSE) {
    res <- bfcquery(bfc, query = hash, field = "rname")
    as.logical(nrow(res))
}

#' @param bfc Object of class BiocFileCache, created by a call to 
#' BiocFileCache::BiocFileCache()
#' @param hash unique hash representing a query.
#' 
#' This function checks if a cache entry is a valid RDS file.
#' Returns TRUE if the cache entry is valid, FALSE otherwise.
#' In the case of an invalid file the cache entry and file are 
#' deleted.
#' @importFrom BiocFileCache bfcremove
#' @keywords Internal
.checkValidCache <- function(bfc, hash) {
    res <- bfcquery(bfc, query = hash, field = "rname")
    if(nrow(res) == 0) {
        return(FALSE)
    } else {
        ## check this is a valid RDS file
        ## remove the cache entry if it's not a valid RDS
        test <- tryCatch(is.list(infoRDS(res$rpath[1])), 
                         error = function(e) { return(FALSE) })
        if(!test) 
            bfcremove(bfc, res$rid[1])
        return(test)
    }
}

.biomartCacheLocation <- function() {
    Sys.getenv(x = "BIOMART_CACHE", 
               unset = rappdirs::user_cache_dir(appname="biomaRt"))
}

biomartCacheClear <- function() {
    cache <- .biomartCacheLocation()
    bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
    removebfc(bfc, ask = FALSE)
}

biomartCacheInfo <- function() {
    cache <- .biomartCacheLocation()
    
    if(!file.exists(cache)) {
        message("biomaRt cache uninitialized\n", 
                "- Location: ", cache)
    } else {
        
        bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
        files <- bfcinfo(bfc)$rpath
        total_size <- sum(file.size(files))
        size_obj <- structure(total_size, class = "object_size")
    
        message("biomaRt cache\n", 
                "- Location: ", cache, "\n",
                "- No. of files: ", length(files), "\n",
                "- Total size: ", format(size_obj, units = "auto"))
    }
    return(invisible(cache))
}

