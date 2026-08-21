.getEnsemblSSL <- function() {
  cache <- .biomartCacheLocation()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  if (.checkInCache(bfc, hash = "ensembl-ssl-settings-httr2")) {
    settings <- .readFromCache(bfc, "ensembl-ssl-settings-httr2")
    if (length(settings) > 0L) {
      .Deprecated(
        package = "biomaRt",
        msg = paste(
          "Modifying the Ensembl SSL settings is deprecated and will be removed in the next release. ",
          "Please get in touch as soon as possible if your workflow requires this function to work."
        )
      )
    }
    return(settings)
  }
  return(list())
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
#' found in [curl::curl_options()].
#' @author Mike Smith
#'
#' @examples
#' \dontrun{
#' ssl_settings <- list(
#'   "ssl_cipher_list" = "DEFAULT@SECLEVEL=1",
#'   "ssl_verifypeer"  = FALSE
#' )
#' setEnsemblSSL(ssl_settings)
#' }
#'
#' @importFrom utils modifyList
#' @export
setEnsemblSSL <- function(settings) {
  .Deprecated(
    package = "biomaRt",
    msg = paste(
      "setEnsemblSSL() is deprecated and will be removed in the next release. ",
      "No modification of SSL settings should be required. ",
      "Please get in touch as soon as possible if your workflow requires this function to work."
    )
  )
  stopifnot(is.list(settings))

  cache <- .biomartCacheLocation()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)

  existing_config <- .getEnsemblSSL()

  if (length(settings) == 0L) {
    updated_config <- list()
  } else {
    updated_config <- modifyList(existing_config, settings)
  }

  .addToCache(
    bfc,
    updated_config,
    hash = "ensembl-ssl-settings-httr2",
    update = TRUE
  )
  return(invisible(TRUE))
}
