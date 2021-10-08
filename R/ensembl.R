## location of Ensembl specific functions

.getArchiveList <- function(https = TRUE, httr_config = list()) {
  
  url_worked <- FALSE
  mirrors <- c("www", "asia", "uswest", "useast")
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

listEnsembl <- function(mart = NULL, version = NULL, 
                        GRCh = NULL, mirror = NULL, verbose = FALSE){
  
  host <- .constructEnsemblURL(mirror = mirror, version = version, GRCh = GRCh)
  port <- ifelse(grepl("https", host)[1], yes = 443, no = 80)
  ensemblRedirect <- is.null(mirror)
  
  httr_config <- .getEnsemblSSL()
  
  marts <- .listMarts(mart = mart, host = host, verbose = verbose, httr_config = httr_config,
                      port = port, ensemblRedirect = ensemblRedirect)
  
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
  
  if(!is.null(mirror) & (!is.null(version) | !is.null(GRCh))){
    warning("version or GRCh arguments cannot be used together with the mirror argument.\n", 
            "We will ignore the mirror argument and connect to the main Ensembl site.",
            call. = FALSE) 
    mirror <- NULL
  }
  
  if(!is.null(version) & !is.null(GRCh)) {
    stop("version or GRCh arguments cannot be used together.\n", 
         "Please specify only the 'version' or 'GRCh' argument.",
         call. = FALSE) 
  }
  
  if(!is.null(version)) {
    archives <- listEnsemblArchives(https = FALSE)
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
    if(!(mirror %in% c("www", "uswest", "useast", "asia"))) {
      warning("Invalid mirror. Select a mirror from [www, uswest, useast, asia].\n",
              "Default when no mirror is specified is to use ",
              "www.ensembl.org which may be automatically redirected." )
      host <- "https://www.ensembl.org"
    } else {
      host <- paste0("https://", mirror, ".ensembl.org")
    }
  }
  
  if(is.null(host)) {
    host = "https://www.ensembl.org"
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

  if(tolower(biomart) == "ensembl" | tolower(biomart) == "genes") {
    biomart = "ENSEMBL_MART_ENSEMBL"
  }
  if(tolower(biomart) == "snp" | tolower(biomart) == "snps"){
    biomart = "ENSEMBL_MART_SNP"
  }
  if(tolower(biomart) == "regulation"){
    biomart = "ENSEMBL_MART_FUNCGEN"
  }
  if(tolower(biomart) == "vega"){
    biomart = "ENSEMBL_MART_VEGA"
  }
  if(tolower(biomart) == "mouse_strains"){
    biomart = "ENSEMBL_MART_MOUSE"
  }
  
  ## create the host URL & turn off redirection if a mirror is specified
  if(missing(host)) {
    host <- .constructEnsemblURL(version = version, GRCh = GRCh, mirror = mirror)
    ensemblRedirect <- is.null(mirror)
  } else {
    ensemblRedirect <- FALSE
  }
  
  ## choose the port based on whether we use https or not
  port <- ifelse(grepl(pattern = "https://", x = host), 
                 yes = 443, no = 80)
  
  ## test https connection and store required settings
  httr_config <- .getEnsemblSSL()
  
  ens <- .useMart(biomart = biomart, 
                  dataset = dataset, 
                  host = host, 
                  verbose = verbose,
                  port = port,
                  ensemblRedirect = ensemblRedirect,
                  httr_config = httr_config)	   
  return(ens)
}

##############################################

listEnsemblGenomes <- function(includeHosts = FALSE){
  
  hosts <- c("https://protists.ensembl.org/",
             "https://fungi.ensembl.org/",
             "https://metazoa.ensembl.org/",
             "https://plants.ensembl.org/")
  
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

useEnsemblGenomes <- function(biomart, dataset) {
  
  if(missing(biomart)) {
    stop("You must provide the argument 'biomart'\n",
         "Available Ensembl Genomes Marts can be viewed with ",
         "the function listEnsemblGenomes()")
  }
  
  marts <- listEnsemblGenomes(includeHosts = TRUE)
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