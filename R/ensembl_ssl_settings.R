.test_ensembl <- function(config = list()) {
  main <- request("https://www.ensembl.org/index.html?redirect=no") |>
    req_timeout(5) |>
    req_options(!!!config)
  useast <- request("https://useast.ensembl.org/index.html?redirect=no") |>
    req_timeout(5) |>
    req_options(!!!config) |>
    ## hopefully temporary work around for 403 error on this mirror
    req_user_agent(
      "Mozilla/5.0 (X11; Linux i686; rv:127.0) Gecko/20100101 Firefox/127.0"
    )

  main |> req_perform()
  useast |> req_perform()
}

.checkEnsemblSSL <- function() {
  ensembl_config <- list()
  test <- try(.test_ensembl(config = ensembl_config), silent = TRUE)
  while (is(test, "try-error")) {
    if (
      grepl(
        test[1], ## This address problems with Ubuntu 20.04 et al and the Ensembl https certificates
        pattern = "sslv3 alert handshake failure",
        fixed = TRUE
      )
    ) {
      ensembl_config[["ssl_cipher_list"]] <- "DEFAULT@SECLEVEL=1"
    } else if (
      grepl(
        x = test[1], ## two reported error messages solved with the same setting
        pattern = "(unable to get local issuer certificate)|(server certificate verification failed)"
      )
    ) {
      ensembl_config[["ssl_verifypeer"]] <- FALSE
    } else if (
      grepl(
        x = test[1], ## We end up here if the test timed out
        pattern = "Timeout was reached",
        fixed = TRUE
      )
    ) {
      ## Time out is unfortunate, but lets not inform the user since it might not be a problem.
      ## Quit the testing and proceed
      break
    } else {
      message(
        "Possible SSL connectivity problems detected.\n",
        "Please report this issue at https://github.com/grimbough/biomaRt/issues\n",
        test[1]
      )
      ## We can't fix this, so just quit the tests
      break
    }
    test <- try(.test_ensembl(config = ensembl_config), silent = TRUE)
  }
  return(ensembl_config)
}

.getEnsemblSSL <- function() {
  cache <- .biomartCacheLocation()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  if (.checkInCache(bfc, hash = "ensembl-ssl-settings-httr2")) {
    ensembl_config <- .readFromCache(bfc, "ensembl-ssl-settings-httr2")
  } else {
    ensembl_config <- .checkEnsemblSSL()
    .addToCache(bfc, ensembl_config, hash = "ensembl-ssl-settings-httr2")
  }
  return(ensembl_config)
}



#' Save system specific SSL settings for contacting Ensembl
#' 
#' On some systems specific SSL settings have to be applied to allow https
#' connections to the Ensembl servers.  This function allows these to be saved
#' in the biomaRt cache, so they will be retrieved each time they are needed.
#' biomaRt will try to determine them automatically, but this function can be
#' used to set them manually if required.
#' 
#' 
#' @param settings A named list. Each entry should be a valid curl option, as
#' found in \code{\link[curl]{curl_options}}.
#' @author Mike Smith
#' @examples
#' 
#' 
#' \dontrun{
#'   ssl_settings <- list("ssl_cipher_list" = "DEFAULT@SECLEVEL=1",
#'                        "ssl_verifypeer"  = FALSE)
#'   setEnsemblSSL(ssl_settings)
#' }
#' 
#' @export setEnsemblSSL
setEnsemblSSL <- function(settings) {
  stopifnot(is.list(settings))

  cache <- .biomartCacheLocation()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)

  existing_config <- .getEnsemblSSL()
  updated_config <- existing_config

  if (length(settings) == 0L) {
    updated_config <- list()
  } else {
    for (i in seq_along(settings)) {
      updated_config[[names(settings)[i]]] <- settings[[i]]
    }
  }

  .addToCache(
    bfc,
    updated_config,
    hash = "ensembl-ssl-settings-httr2",
    update = TRUE
  )
  return(invisible(TRUE))
}
