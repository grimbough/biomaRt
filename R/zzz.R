.onLoad <- function(libname, pkgname) {

    ## This address problems with Ubuntu 20.04 and the Ensembl https certificates
    test <- try(httr::GET("https://www.ensembl.org"), silent = TRUE)
    if(inherits(test, "try-error")) {
      if(grepl(test[1], pattern = "sslv3 alert handshake failure")) {
        httr_config <- httr::config(ssl_verifypeer = FALSE, 
                                    ssl_cipher_list = "DEFAULT@SECLEVEL=1")
        httr::set_config(httr_config)
      }
    }
}
