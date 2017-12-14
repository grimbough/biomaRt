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
    
    tab <- do.call("rbind", extracted)
    colnames(tab) <- c("url", "version", "date")
    tab <- tab[,c(2,3,1)]
    
    return(tab)
}

