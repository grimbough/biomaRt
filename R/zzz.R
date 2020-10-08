.onLoad <- function(libname, pkgname) {
    
    ## This address problems with Ubuntu 20.04 and the Ensembl https certificates
    if( Sys.info()["sysname"] == "Linux" ) {
        httr_config <- httr::config(ssl_verifypeer = FALSE, ssl_cipher_list = "DEFAULT@SECLEVEL=1")
        httr::set_config(httr_config)
    }
    
}
