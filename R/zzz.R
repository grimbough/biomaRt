.test_ensembl <- function() {
  httr::GET("https://www.ensembl.org/index.html", timeout(3), set_cookies(redirect_mirror = "no"))
  httr::GET("https://uswest.ensembl.org/index.html", timeout(3), set_cookies(redirect_mirror = "no"))
}

.check_ensembl_ssl <- function() {
  
  test <- try(stop(), silent = TRUE)
  while(is(test, "try-error")) {
    test <- try(.test_ensembl(), silent = TRUE)
    
    if(is(test, "try-error")) {
      
      if(grepl(test[1], ## This address problems with Ubuntu 20.04 et al and the Ensembl https certificates
               pattern = "sslv3 alert handshake failure")) {
        message("Failed test 1: ", test[1])
        new_config <- httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1")
        
      } else if(grepl(x = test[1], ## two reported error messages solved with the same setting
                      pattern = "(unable to get local issuer certificate)|(server certificate verification failed)")) {
        message("Failed test 2: ", test[1])
        new_config <- httr::config(ssl_verifypeer = FALSE)
      } else {
        message("Unknown error encountered: ", test[1])
        ## We can't fix this, so just quit
        break;
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
