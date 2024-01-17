.test_ensembl <- function(config = list()) {
  main <- request("https://www.ensembl.org/index.html?redirect=no") |>
    req_timeout(5) |>
    req_options(!!!config)
  useast <- request("https://useast.ensembl.org/index.html?redirect=no") |>
    req_timeout(5) |>
    req_options(!!!config)
  
  main |> req_perform()
  useast |> req_perform()
}

.checkEnsemblSSL <- function() {
  
  ensembl_config <- list()
  i <- 1
  test <- try(.test_ensembl(config = ensembl_config), silent = TRUE)
  while(is(test, "try-error")) {
    
    if(grepl(test[1], ## This address problems with Ubuntu 20.04 et al and the Ensembl https certificates
             pattern = "sslv3 alert handshake failure")) {
      ensembl_config[[ "ssl_cipher_list" ]] <- "DEFAULT@SECLEVEL=1"
    } else if (grepl(x = test[1], ## two reported error messages solved with the same setting
                     pattern = "(unable to get local issuer certificate)|(server certificate verification failed)")) {
      ensembl_config[[ "ssl_verifypeer" ]] <- FALSE
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
    test <- try(.test_ensembl(config = ensembl_config), silent = TRUE)
  } 
  return(ensembl_config)
}

.getEnsemblSSL <- function() {
  
  cache <- .biomartCacheLocation()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  if(.checkInCache(bfc, hash = "ensembl-ssl-settings-httr2")) {
    ensembl_config <- .readFromCache(bfc, "ensembl-ssl-settings-httr2")
  } else {
    ensembl_config <- .checkEnsemblSSL()
    .addToCache(bfc, ensembl_config, hash = "ensembl-ssl-settings-httr2")
  }
  return(ensembl_config)
}

setEnsemblSSL <- function(settings) {
  
  stopifnot(is.list(settings))
  
  cache <- .biomartCacheLocation()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  
  existing_config <- .getEnsemblSSL()
  updated_config <- existing_config
  
  if(length(settings) == 0L) {
    updated_config <- list()
  } else {
    for(i in seq_along(settings)) {
      updated_config[[ names(settings)[i] ]] <- settings[[i]]
    }
  }
  
  .addToCache(bfc, updated_config, hash = "ensembl-ssl-settings-httr2", update = TRUE)
  return(invisible(TRUE))
}
