## location of Ensembl specific functions

## scrapes the ensembl website for the list of current archives and returns
## a data frame containing the versions and their URL
listEnsemblArchives <- function() {
    
    html <- htmlParse("http://www.ensembl.org/info/website/archives/index.html?redirect=no")
    
    archive_box <- getNodeSet(html, path = "//div[@class='plain-box float-right archive-box']")[[1]]
    
    archive_box_string <- toString.XMLNode(archive_box)
    
    archives <- strsplit(archive_box_string, split = "<li>")[[1]][-1]
    
    extracted <- str_extract_all(string = archives, 
                    pattern = "Ensembl [A-Za-z0-9 ]{2,6}|http://.*ensembl\\.org|[A-Z][a-z]{2} [0-9]{4}")
    
    ## split the version number into a separate column
    extracted <- lapply(extracted, FUN = function(x) {
        version <- str_match(x[2], pattern = ".+ ([a-zA-Z0-9]+)$")[2]
        return( c(x, version) )
    })
    
    tab <- do.call("rbind", extracted)
    colnames(tab) <- c("url", "name", "date", "version")
    tab <- tab[,c(2,3,1,4)]
    tab[,'url'] <- tolower(tab[,'url'])
    
    return(tab)
}


listEnsembl <- function(mart = NULL, host="www.ensembl.org",version = NULL, GRCh = NULL, mirror = NULL,verbose = FALSE){
    
    if(!is.null(mirror) & (!is.null(version) | !is.null(GRCh))){
        warning("version or GRCh arguments can not be used together with the mirror argument.  Will ignore the mirror argument and connect to default ensembl host") 
        mirror = NULL
    }
    
    if(!is.null(version)){
        host = paste("e",version,".ensembl.org",sep="")
    }
    if(!is.null(GRCh)){
        if(GRCh == 37){ 
            host = paste("grch",GRCh,".ensembl.org",sep="")	
        }
        else{
            print("Only 37 can be specified for GRCh version")
        }
    }
    
    if(!is.null(mirror)){
        if(!(mirror %in% c("www", "uswest", "useast", "asia"))) {
            warning("Invalid mirror select a mirror from [www, uswest, useast, asia].\n",
                    "default when no mirror is specified is to be redirected to ",
                    "www.ensembl.org")
        } else {
            host <- paste0(mirror, ".ensembl.org")
        }
    }
    
    marts = listMarts(mart = mart, host = host, verbose = verbose)
    sel = which(marts$biomart == "ENSEMBL_MART_ENSEMBL")
    if(length(sel) > 0){ 
        marts$biomart[sel] = "ensembl"
    }
    sel = which(marts$biomart == "ENSEMBL_MART_SNP")
    if(length(sel) > 0){ 
        marts$biomart[sel] = "snp"
    }
    sel = which(marts$biomart == "ENSEMBL_MART_FUNCGEN")
    if(length(sel) > 0){ 
        marts$biomart[sel] = "regulation"
    }
    sel = which(marts$biomart == "ENSEMBL_MART_VEGA")
    if(length(sel) > 0){ 
        marts$biomart[sel] = "vega"
    }
    return(marts)
}


useEnsembl <- function(biomart, dataset, host = "www.ensembl.org", version = NULL, GRCh = NULL, mirror = NULL, verbose = FALSE){
    
    if(!is.null(mirror) & (!is.null(version) | !is.null(GRCh))){
        warning("version or GRCh arguments can not be used together with the mirror argument.  Will ignore the mirror argument and connect to default ensembl host") 
        mirror = NULL
    }
    
    if(!is.null(version)){
        archives <- listEnsemblArchives()
        idx <- match(version, archives[,'version'], nomatch = NA)
        if(is.na(idx)) {
            stop('Specified Ensembl version is not available.\n',
                 'Use listEnsemblArchives() to view available versions.',
                 call. = FALSE)
        }
        host = archives[idx, 'url']
    }	   
    
    if(!is.null(GRCh)){
        if(GRCh == 37){
            host = paste("grch",GRCh,".ensembl.org",sep="")
        }
        else{
            print("Only 37 can be specified for GRCh version")
        }
    }
    
    if(!is.null(mirror)){
        if(!(mirror %in% c("www", "uswest", "useast", "asia"))) {
            warning("Invalid mirror select a mirror from [www, uswest, useast, asia].\n",
                    "default when no mirror is specified is to be redirected to ",
                    "www.ensembl.org")
        } else {
            host <- paste0(mirror, ".ensembl.org")
        }
    }
    
    if(biomart == "ensembl"){
        biomart = "ENSEMBL_MART_ENSEMBL"
    }
    if(biomart == "snp"){
        biomart = "ENSEMBL_MART_SNP"
    }
    if(biomart == "regulation"){
        biomart = "ENSEMBL_MART_FUNCGEN"
    }
    if(biomart == "vega"){
        biomart = "ENSEMBL_MART_VEGA"
    }
    
    ens = useMart(biomart = biomart, 
                  dataset = dataset, 
                  host = host, 
                  verbose = verbose)	   
    return(ens)
}
