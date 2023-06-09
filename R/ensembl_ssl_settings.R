.test_ensembl <- function() {
  httr::GET("https://www.ensembl.org/index.html", timeout(5), set_cookies(redirect_mirror = "no"))
  httr::GET("https://useast.ensembl.org/index.html", timeout(5), set_cookies(redirect_mirror = "no"))
}

.checkEnsemblSSL <- function() {
  
  ## we'll modify the global httr setting in this function
  ## get the current settings so we can reset it later
  prev_config <- getOption("httr_config")
  on.exit(options("httr_config" = prev_config))
  
  ensembl_config <- list()
  i <- 1
  test <- try(.test_ensembl(), silent = TRUE)
  while(is(test, "try-error")) {
    
    if(grepl(test[1], ## This address problems with Ubuntu 20.04 et al and the Ensembl https certificates
             pattern = "sslv3 alert handshake failure")) {
      ensembl_config[[i]] <- httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1")
    } else if (grepl(x = test[1], ## two reported error messages solved with the same setting
                     pattern = "(unable to get local issuer certificate)|(server certificate verification failed)")) {
      ensembl_config[[i]] <- new_config <- httr::config(ssl_verifypeer = FALSE)
    } else if (grepl(x = test[1], ## We end up here if the test timed out
                     pattern = "Timeout was reached")) {
      ## Time out is unfortunate, but lets not inform the user since it might not be a problem.
      ## Quit the testing and proceed
      break; 
    } else {
      message("Possible SSL connectivity problems detected.\n",
              "Please report this issue at https://github.com/grimbough/biomaRt/issues\n",
              test[1])
      ## We can't fix this, so just quit the tests
      break;
    }
    httr::set_config(ensembl_config[[i]], override = FALSE)
    i <- i + 1
    test <- try(.test_ensembl(), silent = TRUE)
  } 
  return(ensembl_config)
}

.getEnsemblSSL <- function() {
  
  cache <- .biomartCacheLocation()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  if(.checkInCache(bfc, hash = "ensembl-ssl-settings")) {
    ensembl_config <- .readFromCache(bfc, "ensembl-ssl-settings")
  } else {
    ensembl_config <- .checkEnsemblSSL()
    .addToCache(bfc, ensembl_config, hash = "ensembl-ssl-settings")
  }
  return(ensembl_config)
}