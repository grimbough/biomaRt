## location of Ensembl specific functions

## scrapes the ensembl website for the list of current archives and returns
## a data frame containing the versions and their URL
listEnsemblArchives <- function(https = TRUE) {
  
  url <- ifelse(https,
                "https://www.ensembl.org/info/website/archives/index.html",
                "http://www.ensembl.org/info/website/archives/index.html")
  
  html <- xml2::read_html(GET(url))
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

listEnsembl <- function(mart = NULL, host="www.ensembl.org", version = NULL, GRCh = NULL, mirror = NULL, verbose = FALSE){
  
  host <- .constructEnsemblURL(mirror = mirror, version = version, GRCh = GRCh)
  ensemblRedirect <- is.null(mirror)
  
  marts <- .listMarts(mart = mart, host = host, verbose = verbose, 
                      ensemblRedirect = ensemblRedirect)
  
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
    archives <- listEnsemblArchives()
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
  host <- .constructEnsemblURL(version = version, GRCh = GRCh, mirror = mirror)
  ensemblRedirect = is.null(mirror)
  
  ## choose the port based on whether we use https or not
  port <- ifelse(grepl(pattern = "https://", x = host), 
                 yes = 443, no = 80)
  
  ens = .useMart(biomart = biomart, 
                 dataset = dataset, 
                 host = host, 
                 verbose = verbose,
                 port = port, 
                 ensemblRedirect = ensemblRedirect)	   
  return(ens)
}

##############################################

listEnsemblGenomes <- function(includeHosts = FALSE){
  
  hosts <- c("https://protists.ensembl.org/",
             "https://fungi.ensembl.org/",
             "https://metazoa.ensembl.org/",
             "https://plants.ensembl.org/")
  
  marts <- lapply(hosts, FUN = function(x) { 
    as.data.frame(
      .listMarts(host = x, mart = NULL, 
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
  
  ens <- .useMart(biomart = biomart, 
                 dataset = dataset, 
                 host = host, 
                 verbose = FALSE,
                 port = 443, 
                 ensemblRedirect = FALSE)	  
  
  return(ens)
}