## location of Ensembl specific functions

.getArchiveList <- function(https = TRUE, httr_config = list()) {
  
  url_worked <- FALSE
  mirrors <- c("www", "asia", "useast")
  protocol <- ifelse(https, "https://", "http://" )
  
  while(!url_worked) {
    
    if(length(mirrors) == 0) {
      stop("Unable to contact any Ensembl mirror")
    }

    mirror_option <- mirrors[1]
    mirrors <- mirrors[-1]

    url <- paste0(protocol, mirror_option, ".ensembl.org/info/website/archives/index.html?redirect=no")
  
    html <- GET(url, config = httr_config)
  
    ## this is TRUE if there's an HTTP error or we get the Ensembl error page
    if(identical(status_code(html), 200L) && 
       !grepl("The Ensembl web service you requested is temporarily unavailable", content(html))) {
      return( content(html) )
    }
      
  } 

}

.currentEnsemblVersion <- function() {
    archives <- listEnsemblArchives()
    current <- archives[archives$current_release == "*", ]
    return(current)
}

## scrapes the ensembl website for the list of current archives and returns
## a data frame containing the versions and their URL
listEnsemblArchives <- function(https) {
    
  if(!missing(https)) {
    warning("Ensembl will soon enforce the use of https.\n",
    "As such the 'https' argument will be deprecated in the next release.")
  }
  https <- TRUE
    
  .listEnsemblArchives(https = https, httr_config = list())
}

.listEnsemblArchives <- function(https = TRUE, httr_config) {
  
  html <- .getArchiveList(https, httr_config)
  html <- htmlParse( html )
  
  archive_box <- getNodeSet(html, path = "//div[@class='plain-box float-right archive-box']")[[1]]
  
  archive_box_string <- toString.XMLNode(archive_box)
  
  archives <- strsplit(archive_box_string, split = "<li>")[[1]][-1]
  
  extracted <- str_extract_all(string = archives, 
                               pattern = "Ensembl [A-Za-z0-9 ]{2,6}|http[s]?://.*ensembl\\.org|[A-Z][a-z]{2} [0-9]{4}")
  
  ## split the version number into a separate column
  extracted <- lapply(extracted, FUN = function(x) {
    version <- str_match(x[2], pattern = ".+ ([a-zA-Z0-9]+)$")[2]
    return( c(x, version) )
  })
  
  current <- ifelse(stringr::str_detect(archives, "- this site"), "*", "")
  
  tab <- do.call("rbind", extracted)
  tab <- cbind(tab, current)
  
  dframe <- data.frame("name" = as.character(tab[,2]), 
                       "date" = as.character(tab[,3]), 
                       "url" = stringr::str_replace(tolower(as.character(tab[,1])),
                                                    "http://",
                                                    "https://"),
                       "version" = as.character(tab[,4]), 
                       "current_release" = as.character(tab[,5]),
                       stringsAsFactors = FALSE)
  return(dframe)
}

.listEnsembl <- function(mart = NULL, version = NULL, GRCh = NULL, 
                         mirror = NULL, verbose = FALSE) {
    
    if(is.null(version)) {
        version_num <- .currentEnsemblVersion()["version"]
    } else {
        version_num <- version
    }
    
    ## determine if a cached version exists and if it's less than one week old
    cache <- .biomartCacheLocation()
    bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
    use_cached_version <- FALSE
    if(.checkInCache(bfc, hash = paste0("ensembl-marts-", version_num))) {
        cache_entry <- bfcquery(x = bfc, query = paste0("ensembl-marts-", version_num))
        if( (nrow(cache_entry) == 1) && (as.Date(Sys.time()) - as.Date(cache_entry$create_time) < 7) ) {
            use_cached_version <- TRUE
        } else {
            bfcremove(bfc, cache_entry$rid)
        }
    }
    
    if(use_cached_version) {
        marts <- .readFromCache(bfc, paste0("ensembl-marts-", version_num))
    } else {
        host <- .constructEnsemblURL(mirror = mirror, version = version, GRCh = GRCh)
        port <- ifelse(grepl("https", host)[1], yes = 443, no = 80)
        ensemblRedirect <- is.null(mirror)
        
        httr_config <- .getEnsemblSSL()
        
        marts <- .listMarts(mart = mart, host = host, verbose = verbose, httr_config = httr_config,
                            port = port, ensemblRedirect = ensemblRedirect)
        
        .addToCache(bfc, marts, hash = paste0("ensembl-marts-", version_num))
    }
    
    return(marts)
    
}
    

listEnsembl <- function(mart = NULL, version = NULL, GRCh = NULL, 
                        mirror = NULL, verbose = FALSE) {
    
  marts <- .listEnsembl(mart = mart, version = version, GRCh = GRCh,
                        mirror = mirror, verbose = verbose)
  
  sel = which(marts$biomart == "ENSEMBL_MART_ENSEMBL")
  if(length(sel) > 0){ 
    marts$biomart[sel] = "genes"
  }
  sel = which(marts$biomart == "ENSEMBL_MART_SNP")
  if(length(sel) > 0){ 
    marts$biomart[sel] = "snps"
  }
  sel = which(marts$biomart == "ENSEMBL_MART_FUNCGEN")
  if(length(sel) > 0){ 
    marts$biomart[sel] = "regulation"
  }
  sel = which(marts$biomart == "ENSEMBL_MART_VEGA")
  if(length(sel) > 0){ 
    marts$biomart[sel] = "vega"
  }
  sel = which(marts$biomart == "ENSEMBL_MART_MOUSE")
  if(length(sel) > 0){ 
    marts$biomart[sel] = "mouse_strains"
  }
  return(marts)
}

#' creates an Ensembl URL based on the arguments provided to useEnsembl.
#' If there are conflicting options, order of precedence is:
#' GRCh, version, mirror
#' Default return value is https://www.ensembl.org
.constructEnsemblURL <- function(mirror = NULL, version = NULL, GRCh = NULL) {
  
  host <- NULL
  
  if(!is.null(mirror) && (!is.null(version) || !is.null(GRCh))){
    warning("version or GRCh arguments cannot be used together with the mirror argument.\n", 
            "We will ignore the mirror argument and connect to the main Ensembl site.",
            call. = FALSE) 
    mirror <- NULL
  }
  
  if(!is.null(version) && !is.null(GRCh)) {
    stop("version or GRCh arguments cannot be used together.\n", 
         "Please specify only the 'version' or 'GRCh' argument.",
         call. = FALSE) 
  }
  
  if(!is.null(version)) {
    archives <- .listEnsemblArchives(https = TRUE, httr_config = list())
    idx <- match(version, archives[,'version'], nomatch = NA)
    if(is.na(idx)) {
      stop("Specified Ensembl version is not available.\n",
           "Use listEnsemblArchives() to view available versions.",
           call. = FALSE)
    }
    host <- archives[idx, 'url']
  }	  
  
  if(!is.null(GRCh)){
    if(GRCh == 37){
      host <- paste0("https://grch", GRCh, ".ensembl.org")
    } else {
      warning("Only 37 can be specified for GRCh version. Using the current version.",
              call. = FALSE)
    }
  }
  
  if(!is.null(mirror)){
    if(!(mirror %in% c("www", "useast", "asia"))) {
      warning("Invalid mirror. Select a mirror from [www, useast, asia].\n",
              "Default when no mirror is specified is to use ",
              "www.ensembl.org which may be automatically redirected." )
      host <- "https://www.ensembl.org"
    } else {
      host <- paste0("https://", mirror, ".ensembl.org")
    }
  }
  
  if(is.null(host)) {
    host <- "https://www.ensembl.org"
  }
  
  return(host)
  
}

useEnsembl <- function(biomart, dataset, host, 
                       version = NULL, GRCh = NULL, mirror = NULL, verbose = FALSE){
  
  if(missing(biomart)) {
    stop("You must provide the argument 'biomart'\n",
         "Available Ensembl Marts can be viewed with ",
         "the function listEnsembl()")
  }

  biomart <- switch (tolower(biomart),
      "ensembl" = "ENSEMBL_MART_ENSEMBL",
      "genes"   = "ENSEMBL_MART_ENSEMBL",
      "snp"     = "ENSEMBL_MART_SNP",
      "snps"    = "ENSEMBL_MART_SNP",
      "regulation" = "ENSEMBL_MART_FUNCGEN",
      "mouse_strains" = "ENSEMBL_MART_MOUSE",
      "vega"          = "ENSEMBL_MART_VEGA",
      biomart
  )
  
  ## test https connection and store required settings
  httr_config <- .getEnsemblSSL()
  
  ## a crude check to ensure the sub-domain is included.  Otherwise queries will fail
  if(!missing(host)) {
    no_subdomain <- grepl(x = host, pattern = "http[s]?://ensembl", fixed = FALSE)
  } else {
    no_subdomain <- FALSE
  }
  
  ## create the host URL & turn off redirection if a mirror is specified
  if(missing(host) || no_subdomain ) {
    
    if(no_subdomain) {
      warning("You cannot use the host 'ensembl.org'.\n",
              "Please provide a subdomain e.g. www.ensembl.org or use one of the 'mirror', 'version', 'GRCh' arguments")
    }
    
    if(is.null(version) && is.null(GRCh)) {  
        mirror <- .chooseEnsemblMirror(mirror = mirror, httr_config = httr_config)
    }
    host <- .constructEnsemblURL(version = version, GRCh = GRCh, mirror = mirror)
    ensemblRedirect <- is.null(mirror)
  } else {
    ensemblRedirect <- FALSE
  }
  
  ## choose the port based on whether we use https or not
  port <- ifelse(grepl(pattern = "https://", x = host), 
                 yes = 443, no = 80)
  
  if(grepl(x = host, pattern = "www|useast|asia")) {
    marts <- .listEnsembl(version = version, GRCh = GRCh, mirror = mirror)
  } else {
    marts <- .listMarts(host = host, port = port, httr_config = httr_config, ensemblRedirect = FALSE)
  }
  
  mindex = NA
  if(!missing(biomart)){ 
      mindex=match(biomart,marts$biomart)
  }
  if(is.na(mindex))
      stop("Incorrect BioMart name, use the listMarts function to see which BioMart databases are available")
  
  ## adding option to force use of specified host with ensembl
  redirect <- ifelse(!ensemblRedirect && grepl(x = host, pattern = "ensembl.org"), 
                     "?redirect=no",
                     "")
  
  mart <- Mart( 
      biomart = biomart,
      vschema = "default", 
      host = paste0(host, ":", 
                    port,
                    "/biomart/martservice",
                    redirect),
      httr_config = httr_config
  )
  
  if(grepl("archive", martHost(mart))) {
      
      ## hack to work around redirection of most recent mirror URL
      archives <- .listEnsemblArchives(https = TRUE, httr_config = httr_config)
      current_release <- archives[archives$current_release == "*", 'url']
      if(grepl(martHost(mart), pattern = current_release)) {
          martHost(mart) <- stringr::str_replace(martHost(mart), pattern = current_release, "https://www.ensembl.org")
          martHost(mart) <- stringr::str_replace(martHost(mart), pattern = ":80/", ":443/")
      }
  }
  
  if(!missing(dataset)){
      mart = useDataset(mart = mart, dataset=dataset, verbose = verbose)
  }
  return(mart)
  
}


##############################################

listEnsemblGenomes <- function(includeHosts = FALSE, host = NULL){
  
  ## use the default websites unless an alternative is provided
  if(is.null(host)) {
    hosts <- c("https://protists.ensembl.org/",
               "https://fungi.ensembl.org/",
               "https://metazoa.ensembl.org/",
               "https://plants.ensembl.org/")
  } else {
    hosts <- host
  }
  
  httr_config <- .getEnsemblSSL()
  
  marts <- lapply(hosts, FUN = function(x) { 
    as.data.frame(
      .listMarts(host = x, mart = NULL, httr_config = httr_config,
                 verbose = FALSE, ensemblRedirect = FALSE, 
                 port = 443, includeHosts = includeHosts)
    ) } 
  )
  
  marts <- do.call("rbind", marts)
  
  return(marts)
}

useEnsemblGenomes <- function(biomart, dataset, host = NULL) {
  
  if(missing(biomart)) {
    stop("You must provide the argument 'biomart'\n",
         "Available Ensembl Genomes Marts can be viewed with ",
         "the function listEnsemblGenomes()")
  }
  
  marts <- listEnsemblGenomes(includeHosts = TRUE, host = host)
  if(!biomart %in% marts$biomart) {
    stop(biomart, " is not in the list of available Marts'\n",
         "Available Ensembl Genomes Marts can be viewed with ",
         "the function listEnsemblGenomes()")
  } else {
    martDetails <- marts[which(marts$biomart == biomart), ]
  }
  
  host <- paste0("https://", martDetails$host)
  
  httr_config <- .getEnsemblSSL()
  
  ens <- .useMart(biomart = biomart, 
                 dataset = dataset, 
                 host = host, 
                 verbose = FALSE,
                 port = 443, 
                 ensemblRedirect = FALSE,
                 httr_config = httr_config)	  
  
  return(ens)
}


#' This function submits a small test query to identify a working Ensembl mirror.
#' If no mirror argument is provided it will use "www" as its first choice.  
#' If the selected mirror returns a success (http 200) response it will be used
#' Otherwise another mirror is selected at random and used instead.
#' If all mirrors fail it will return an error
.chooseEnsemblMirror <- function(mirror, httr_config) {
    
    mirrors <- c("www", "asia", "useast")
    
    if(missing(httr_config)) {
        httr_config <- do.call(c, .getEnsemblSSL())
    }
    if(is.list(httr_config)) {
        httr_config <- do.call(c, httr_config)
    }
    
    example_query <- '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "ensembl_gene_id" value = "ENSG00000000003"/>
		<Attribute name = "ensembl_gene_id" />
	</Dataset>
</Query>'
    
    ## create Ensembl URL and stop any redirection to a mirror
    host <- .constructEnsemblURL(mirror = mirror)
    host <- paste0(host, "/biomart/martservice?redirect=no")
    mirror <- str_match(host, pattern = "://([a-z]{3,6})\\.")[1,2]
    
    result <- tryCatch(httr::POST(url = host,
                                  body = list('query' = example_query),
                                  config = httr_config,
                                  content_type("text/plain"), 
                                  timeout(10)),
                       error = function(c) { "timeout" } )
    
    tryAgain <- any(result == "timeout") || httr::status_code(result) == 500
    
    if(tryAgain) { ## try an alternative mirror if ensembl returns 500
        remaining_mirrors <- setdiff(mirrors, mirror)
        while((length(remaining_mirrors) > 0) && (tryAgain)) {
            mirror <- sample(remaining_mirrors, size = 1)
            message("Ensembl site unresponsive, trying ", mirror, " mirror")
            host <- str_replace(host, 
                                pattern = "://([a-z]{3,6})\\.", 
                                replacement = paste0("://", mirror, "."))
            result <- tryCatch(httr::POST(url = host,
                                          body = list('query' = example_query),
                                          config = httr_config,
                                          content_type("text/plain"), 
                                          timeout(10)),
                               error = function(c) { "timeout" } )
            tryAgain <- any(result == "timeout") || httr::status_code(result) == 500
            if(tryAgain) {
                remaining_mirrors <- setdiff(remaining_mirrors, mirror)
            }
        }
    }
    if(tryAgain) {
        stop("Unable to query any Ensembl site")
    }
    
    return(mirror)
}
