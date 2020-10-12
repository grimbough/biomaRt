.onLoad <- function(libname, pkgname) {

  default_sec <- FALSE
    
    ## This address problems with Ubuntu 20.04 and the Ensembl https certificates
    test <- try(httr::GET("https://www.ensembl.org"), silent = TRUE)
    if(inherits(test, "try-error")) {
      if(grepl(test[1], pattern = "sslv3 alert handshake failure")) {
        default_sec <- TRUE
        new_config <- httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1")
        httr::set_config(new_config, override = TRUE)
      }
    }
  
    ## This address problems with Ubuntu 20.04 and the Ensembl https certificates
    test <- try(httr::GET("https://www.ensembl.org"), silent = TRUE)
    if(inherits(test, "try-error")) {
      if(grepl(test[1], pattern = "unable to get local issuer certificate")) {
        if(default_sec) {
          new_config <- httr::config(ssl_verifypeer = FALSE, 
                                     ssl_cipher_list = "DEFAULT@SECLEVEL=1")
        } else {
          new_config <- httr::config(ssl_verifypeer = FALSE)
        }
        verify_peer <- TRUE
        httr::set_config(new_config, override = FALSE)
      }
    }

}
