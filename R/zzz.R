.test_ensembl <- function() {
  httr::GET("https://www.ensembl.org/index.html", timeout(3), set_cookies(redirect_mirror = "no"))
  httr::GET("https://uswest.ensembl.org/index.html", timeout(3), set_cookies(redirect_mirror = "no"))
}

.check_ensembl_ssl <- function() {
  for(i in seq_len(2)) {
    test <- try(.test_ensembl(), silent = TRUE)
    
    if(inherits(test, "try-error")) {
      
      if(grepl(test[1], ## This address problems with Ubuntu 20.04 et al and the Ensembl https certificates
               pattern = "sslv3 alert handshake failure")) {
        new_config <- httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1")
        
      } else if(grepl(x = test[1], ## two reported error messages solved with the same setting
                      pattern = "(unable to get local issuer certificate)|(server certificate verification failed)")) {
        new_config <- httr::config(ssl_verifypeer = FALSE)
      }
      httr::set_config(new_config, override = FALSE)
    } else {
      ## no need to loop twice if there's no error
      break;
    }
  }
}


.onLoad <- function(libname, pkgname) {
  
  .check_ensembl_ssl()
  
}
