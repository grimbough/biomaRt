.test_ensembl <- function() {
  httr::GET("https://www.ensembl.org/index.html", timeout(3), set_cookies(redirect_mirror = "no"))
  httr::GET("https://uswest.ensembl.org/index.html", timeout(3), set_cookies(redirect_mirror = "no"))
}

.check_ensembl_ssl <- function() {
  
  test <- try(stop(), silent = TRUE)
  test1_done_already <- test2_done_already <- FALSE
  while(is(test, "try-error")) {
    test <- try(.test_ensembl(), silent = TRUE)
    
    if(is(test, "try-error")) {
      
      if(grepl(test[1], ## This address problems with Ubuntu 20.04 et al and the Ensembl https certificates
               pattern = "sslv3 alert handshake failure") && 
        !test1_done_already) {
        new_config <- httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1")
        test1_done_already <- TRUE
      } else if (grepl(x = test[1], ## two reported error messages solved with the same setting
                      pattern = "(unable to get local issuer certificate)|(server certificate verification failed)") &&
                 !test2_done_already) {
        new_config <- httr::config(ssl_verifypeer = FALSE)
      } else if (grepl(x = test[1], ## We end up here if the test timed out
                       pattern = "Timeout was reached")) {
        ## Time out is unfortunate, but lets not inform the user since it might not be a problem.
        ## Quit the testing and proceed
        break; 
      } else {
        message("Possible Ensembl SSL connectivity problems detected.\n",
                "Please see the 'Connection Troubleshooting' section of the biomaRt vignette\n", 
                "vignette('accessing_ensembl', package = 'biomaRt')",
                test[1])
        ## We can't fix this, so just quit the tests
        break
      }
      httr::set_config(new_config, override = FALSE)
    } else {
      ## no need to loop again if there's no error
      break
    }
  }
}


.onLoad <- function(libname, pkgname) {
  .check_ensembl_ssl()
}
