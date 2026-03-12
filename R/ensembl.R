## location of Ensembl specific functions

.checkArchiveList <- function(http_config = list()) {
  ## determine if a cached version exists and if it's less than one week old
  cache <- .biomartCacheLocation()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  cache_entry <- "ensembl-archive-html"
  use_cached_version <- .useCache(
    bfc = bfc,
    cacheEntry = cache_entry,
    numDays = 7L
  )

  if (use_cached_version) {
    archive_html <- .readFromCache(bfc, cache_entry)
  } else {
    archive_html <- .getArchiveList(http_config = http_config)
    .addToCache(bfc, archive_html, hash = cache_entry)
  }

  return(archive_html)
}

#' @importFrom httr2 req_error req_options req_perform req_retry req_timeout request resp_body_string resp_status
.getArchiveList <- function(http_config = list()) {
  mirrors <- c("www", "asia", "useast")

  while (length(mirrors) > 0) {
    url <- paste0(
      "https://",
      mirrors[1],
      ".ensembl.org/info/website/archives/index.html?redirect=no"
    )

    html_request <- request(url) |>
      req_user_agent(
        .biomaRt_user_agent()
      ) |>
      req_timeout(10) |>
      req_options(!!!http_config) |>
      req_retry(max_tries = 3) |>
      req_error(is_error = \(resp) FALSE)

    html <- req_perform(html_request)

    ## this is TRUE if there's an HTTP error or we get the Ensembl error page
    if (
      identical(resp_status(html), 200L) &&
        !grepl(
          "The Ensembl web service you requested is temporarily unavailable",
          resp_body_string(html),
          fixed = TRUE
        )
    ) {
      return(resp_body_string(html))
    }
    mirrors <- mirrors[-1]
  }
  stop("Unable to contact any Ensembl mirror")
}

.currentEnsemblVersion <- function() {
  archives <- listEnsemblArchives()
  current <- archives[archives$current_release == "*", ]
  return(current)
}

## scrapes the ensembl website for the list of current archives and returns
## a data frame containing the versions and their URL

#' Lists the available archived versions of Ensembl
#'
#' Returns a table containing the available archived versions of Ensembl, along
#' with the dates they were created and the URL used to access them.
#'
#' @author Mike Smith
#' @keywords methods
#'
#' @examples
#' listEnsemblArchives()
#'
#' @export
listEnsemblArchives <- function() {
  .listEnsemblArchives(http_config = list())
}

#' @importFrom stringr str_extract_all str_match
.listEnsemblArchives <- function(http_config) {
  html <- .checkArchiveList(http_config)
  html <- xml2::read_html(html)

  archive_box <- as.character(
    xml2::xml_find_first(
      x = html,
      xpath = "//div[contains(@class,'archive-box')]"
    )
  )

  archives <- strsplit(archive_box, split = "<li>", fixed = TRUE)[[1]][-1]

  extracted <- str_extract_all(
    string = archives,
    pattern = "Ensembl [A-Za-z0-9 ]{2,6}|https?://.*ensembl\\.org|[A-Z][a-z]{2} [0-9]{4}"
  )

  ## split the version number into a separate column
  extracted <- lapply(extracted, FUN = function(x) {
    version <- str_match(x[2], pattern = ".+ ([a-zA-Z0-9]+)$")[2]
    return(c(x, version))
  })

  current <- ifelse(
    stringr::str_detect(archives, stringr::fixed("- this site")),
    "*",
    ""
  )

  tab <- do.call("rbind", extracted)
  tab <- cbind(tab, current)

  dframe <- data.frame(
    "name" = as.character(tab[, 2]),
    "date" = as.character(tab[, 3]),
    "url" = stringr::str_replace(
      tolower(as.character(tab[, 1])),
      stringr::fixed("http://"),
      "https://"
    ),
    "version" = as.character(tab[, 4]),
    "current_release" = as.character(tab[, 5]),
    stringsAsFactors = FALSE
  )
  return(dframe)
}

.listEnsembl <- function(
  mart = NULL,
  version = NULL,
  GRCh = NULL,
  mirror = NULL,
  verbose = FALSE
) {
  if (is.null(version)) {
    version_num <- .currentEnsemblVersion()["version"]
  } else {
    version_num <- version
  }

  ## determine if a cached version exists and if it's less than one week old
  cache <- .biomartCacheLocation()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  use_cached_version <- .useCache(
    bfc = bfc,
    cacheEntry = paste0("ensembl-marts-", version_num),
    numDays = 7L
  )

  if (use_cached_version) {
    marts <- .readFromCache(bfc, paste0("ensembl-marts-", version_num))
  } else {
    host <- .constructEnsemblURL(
      mirror = mirror,
      version = version,
      GRCh = GRCh
    )
    port <- ifelse(startsWith(host, "http://")[1], yes = 80, no = 443)
    ensemblRedirect <- is.null(mirror)

    http_config <- .getEnsemblSSL()

    marts <- .listMarts(
      mart = mart,
      host = host,
      verbose = verbose,
      http_config = http_config,
      port = port,
      ensemblRedirect = ensemblRedirect
    )

    .addToCache(bfc, marts, hash = paste0("ensembl-marts-", version_num))
  }

  return(marts)
}


#' lists the available BioMart databases hosted by Ensembl
#'
#' This function returns a list of BioMart databases hosted by Ensembl.  To
#' establish a connection use the [useEnsembl()] function.
#'
#' @param mart mart object created with the useEnsembl function.  This is
#' optional, as you usually use [listMarts()] to see which marts
#' there are to connect to.
#' @param version Ensembl version to connect to when wanting to connect to an
#' archived Ensembl version
#' @param GRCh GRCh version to connect to if not the current GRCh38, currently
#' this can only be 37
#' @param mirror Specify an Ensembl mirror to connect to.  The valid options
#' here are 'www', 'useast', 'asia'.  If no mirror is specified the primary
#' site at www.ensembl.org will be used.
#' @param verbose Give detailed output of what the method is doing, for
#' debugging purposes
#' @param includeHosts If this option is set to `TRUE` a more detailed
#' output is produced, including the URL used to access the corresponding mart.
#' @param host Host to connect to. Use this argument to specify and archive
#' site for [listEnsemblGenomes()] to work with.
#' @author Steffen Durinck, Mike L. Smith
#' @keywords methods
#'
#' @examplesIf interactive()
#' listEnsembl()
#'
#' ## list the default Ensembl Genomes marts
#' listEnsemblGenomes()
#'
#' ## list only the marts available in the Ensmbl Plants 56 archive
#' listEnsemblGenomes(host = "https://eg56-plants.ensembl.org/")
#'
#' @export
listEnsembl <- function(
  mart = NULL,
  version = NULL,
  GRCh = NULL,
  mirror = NULL,
  verbose = FALSE
) {
  marts <- .listEnsembl(
    mart = mart,
    version = version,
    GRCh = GRCh,
    mirror = mirror,
    verbose = verbose
  )

  sel <- which(marts$biomart == "ENSEMBL_MART_ENSEMBL")
  if (length(sel) > 0) {
    marts$biomart[sel] <- "genes"
  }
  sel <- which(marts$biomart == "ENSEMBL_MART_SNP")
  if (length(sel) > 0) {
    marts$biomart[sel] <- "snps"
  }
  sel <- which(marts$biomart == "ENSEMBL_MART_FUNCGEN")
  if (length(sel) > 0) {
    marts$biomart[sel] <- "regulation"
  }
  sel <- which(marts$biomart == "ENSEMBL_MART_VEGA")
  if (length(sel) > 0) {
    marts$biomart[sel] <- "vega"
  }
  sel <- which(marts$biomart == "ENSEMBL_MART_MOUSE")
  if (length(sel) > 0) {
    marts$biomart[sel] <- "mouse_strains"
  }
  return(marts)
}

## creates an Ensembl URL based on the arguments provided to useEnsembl.
## If there are conflicting options, order of precedence is:
## GRCh, version, mirror
## Default return value is https://www.ensembl.org
.constructEnsemblURL <- function(mirror = NULL, version = NULL, GRCh = NULL) {
  host <- NULL

  if (!is.null(mirror) && (!is.null(version) || !is.null(GRCh))) {
    warning(
      "version or GRCh arguments cannot be used together with the mirror argument.\n",
      "We will ignore the mirror argument and connect to the main Ensembl site.",
      call. = FALSE
    )
    mirror <- NULL
  }

  if (!is.null(version) && !is.null(GRCh)) {
    stop(
      "version or GRCh arguments cannot be used together.\n",
      "Please specify only the 'version' or 'GRCh' argument.",
      call. = FALSE
    )
  }

  if (!is.null(version)) {
    archives <- .listEnsemblArchives(http_config = list())
    idx <- match(version, archives[, "version"], nomatch = NA)
    if (is.na(idx)) {
      stop(
        "Specified Ensembl version is not available.\n",
        "Use listEnsemblArchives() to view available versions.",
        call. = FALSE
      )
    }
    host <- archives[idx, "url"]
  }

  if (!is.null(GRCh)) {
    if (GRCh == 37) {
      host <- paste0("https://grch", GRCh, ".ensembl.org")
    } else {
      warning(
        "Only 37 can be specified for GRCh version. Using the current version.",
        call. = FALSE
      )
    }
  }

  if (!is.null(mirror)) {
    if (mirror %in% c("www", "useast", "asia")) {
      host <- paste0("https://", mirror, ".ensembl.org")
    } else {
      warning(
        "Invalid mirror. Select a mirror from [www, useast, asia].\n",
        "Default when no mirror is specified is to use ",
        "www.ensembl.org which may be automatically redirected."
      )
      host <- "https://www.ensembl.org"
    }
  }

  if (is.null(host)) {
    host <- "https://www.ensembl.org"
  }

  return(host)
}


#' Connects to the selected BioMart database and dataset hosted by Ensembl
#'
#' A first step in using the biomaRt package is to select a BioMart database
#' and dataset to use.  The [useEnsembl()] function enables one to connect
#' to a specified BioMart database and dataset hosted by Ensembl without having
#' to specify the Ensembl URL.  To know which BioMart databases are available
#' see the [listEnsembl()] and [listEnsemblGenomes()]
#' functions.  To know which datasets are available within a BioMart database,
#' first select the BioMart database using [useEnsembl()] and then use the
#' [listDatasets()] function on the selected Mart object.
#'
#' The `mirror` argument can be considered as a "preferred choice" when
#' connecting to Ensembl.  If the argument is provided then connectivity to
#' that mirror will be tested.  If it responds positively then the requested
#' mirror will be used.  If the response is a failure each of the remaining
#' mirrors will be selected at random and tested until a working server is
#' found.  Once identified that Ensembl server will be associated with the
#' returned `Mart` object and will be used for all queries.
#'
#' @param biomart BioMart database name you want to connect to. Possible
#' database names can be retrieved with the function [listEnsembl()]
#' @param dataset Dataset you want to use.  To see the different datasets
#' available within a biomaRt you can e.g. do: mart = useEnsembl('genes'),
#' followed by listDatasets(mart).
#' @param host Host to connect to.  Only needs to be specified if different
#' from www.ensembl.org.  For [useEnsemblGenomes()] this argument can be
#' used to specify an archive site.
#' @param version Ensembl version to connect to when wanting to connect to an
#' archived Ensembl version
#' @param GRCh GRCh version to connect to if not the current GRCh38, currently
#' this can only be 37
#' @param mirror Specify an Ensembl mirror to connect to.  The valid options
#' here are 'www', 'useast', 'asia'.  If no mirror is specified the primary
#' site at www.ensembl.org will be used.  Mirrors are not available for the
#' Ensembl Genomes databases.
#' @param verbose Give detailed output of what the method is doing while in
#' use, for debugging
#' @author Steffen Durinck & Mike Smith
#' @keywords methods
#'
#' @examplesIf interactive()
#' mart <- useEnsembl("ENSEMBL_MART_ENSEMBL")
#'
#' ## using the US East mirror
#' us_mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", mirror = "useast")
#'
#' ## using the Arabidopsis thaliana genes dataset in Ensembl Plants
#' plants_mart <- useEnsemblGenomes(
#'   biomart = "plants_mart",
#'   dataset = "athaliana_eg_gene"
#' )
#'
#' ## using the Cucumis melo genes dataset in the Ensembl Plants 56 archive
#' plants_mart <- useEnsemblGenomes(
#'   biomart = "plants_mart",
#'   dataset = "cmelo_eg_gene",
#'   host = "https://feb2023-plants.ensembl.org/"
#' )
#'
#' @export
useEnsembl <- function(
  biomart,
  dataset,
  host,
  version = NULL,
  GRCh = NULL,
  mirror = NULL,
  verbose = FALSE
) {
  if (missing(biomart)) {
    stop(
      "You must provide the argument 'biomart'\n",
      "Available Ensembl Marts can be viewed with ",
      "the function listEnsembl()"
    )
  }

  biomart <- switch(
    tolower(biomart),
    "ensembl" = "ENSEMBL_MART_ENSEMBL",
    "genes" = "ENSEMBL_MART_ENSEMBL",
    "snp" = "ENSEMBL_MART_SNP",
    "snps" = "ENSEMBL_MART_SNP",
    "regulation" = "ENSEMBL_MART_FUNCGEN",
    "mouse_strains" = "ENSEMBL_MART_MOUSE",
    "vega" = "ENSEMBL_MART_VEGA",
    biomart
  )

  ## test https connection and store required settings
  http_config <- .getEnsemblSSL()

  ## a crude check to ensure the sub-domain is included.  Otherwise queries will fail
  if (!missing(host)) {
    no_subdomain <- grepl(
      x = host,
      pattern = "https?://ensembl",
      fixed = FALSE
    )
  } else {
    no_subdomain <- FALSE
  }

  ## create the host URL & turn off redirection if a mirror is specified
  if (missing(host) || no_subdomain) {
    if (no_subdomain) {
      warning(
        "You cannot use the host 'ensembl.org'.\n",
        "Please provide a subdomain e.g. www.ensembl.org or use one of the 'mirror', 'version', 'GRCh' arguments"
      )
    }

    if (is.null(version) && is.null(GRCh)) {
      mirror <- .chooseEnsemblMirror(mirror = mirror, http_config = http_config)
    }
    host <- .constructEnsemblURL(
      version = version,
      GRCh = GRCh,
      mirror = mirror
    )
    ensemblRedirect <- is.null(mirror)
  } else {
    ensemblRedirect <- FALSE
  }

  ## choose the port based on whether we use https or not
  port <- ifelse(startsWith(host, "http://"), yes = 80, no = 443)

  if (grepl(x = host, pattern = "www|useast|asia")) {
    marts <- .listEnsembl(version = version, GRCh = GRCh, mirror = mirror)
  } else {
    marts <- .listMarts(
      host = host,
      port = port,
      http_config = http_config,
      ensemblRedirect = FALSE
    )
  }

  mindex <- NA
  if (!missing(biomart)) {
    mindex <- match(biomart, marts$biomart)
  }
  if (is.na(mindex)) {
    stop(
      "Incorrect BioMart name, use the listMarts function to see which BioMart databases are available"
    )
  }

  ## adding option to force use of specified host with ensembl
  redirect <- ifelse(
    !ensemblRedirect && grepl(x = host, pattern = "ensembl.org"),
    "?redirect=no",
    ""
  )

  mart <- Mart(
    biomart = biomart,
    vschema = "default",
    host = paste0(host, ":", port, "/biomart/martservice", redirect),
    http_config = http_config
  )

  if (grepl("archive", martHost(mart), fixed = TRUE)) {
    ## hack to work around redirection of most recent mirror URL
    archives <- .listEnsemblArchives(http_config = http_config)
    current_release <- archives[archives$current_release == "*", "url"]
    if (grepl(martHost(mart), pattern = current_release)) {
      martHost(mart) <- stringr::str_replace(
        martHost(mart),
        pattern = current_release,
        "https://www.ensembl.org"
      )
      martHost(mart) <- stringr::str_replace(
        martHost(mart),
        pattern = stringr::fixed(":80/"),
        ":443/"
      )
    }
  }

  if (!missing(dataset)) {
    mart <- useDataset(mart = mart, dataset = dataset, verbose = verbose)
  }
  return(mart)
}


##############################################
#' @rdname listEnsembl
#' @export
listEnsemblGenomes <- function(includeHosts = FALSE, host = NULL) {
  ## use the default websites unless an alternative is provided
  if (is.null(host)) {
    host <- c(
      "https://protists.ensembl.org/",
      "https://fungi.ensembl.org/",
      "https://metazoa.ensembl.org/",
      "https://plants.ensembl.org/"
    )
  }

  http_config <- .getEnsemblSSL()

  marts <- lapply(host, FUN = function(x) {
    as.data.frame(
      .listMarts(
        host = x,
        mart = NULL,
        http_config = http_config,
        verbose = FALSE,
        ensemblRedirect = FALSE,
        port = 443,
        includeHosts = includeHosts
      )
    )
  })

  marts <- do.call("rbind", marts)

  return(marts)
}

#' @rdname useEnsembl
#' @export
useEnsemblGenomes <- function(biomart, dataset, host = NULL) {
  if (missing(biomart)) {
    stop(
      "You must provide the argument 'biomart'\n",
      "Available Ensembl Genomes Marts can be viewed with ",
      "the function listEnsemblGenomes()"
    )
  }

  marts <- listEnsemblGenomes(includeHosts = TRUE, host = host)
  if (!biomart %in% marts$biomart) {
    stop(
      biomart,
      " is not in the list of available Marts'\n",
      "Available Ensembl Genomes Marts can be viewed with ",
      "the function listEnsemblGenomes()"
    )
  }
  martDetails <- marts[which(marts$biomart == biomart), ]

  host <- paste0("https://", martDetails$host)

  http_config <- .getEnsemblSSL()

  ens <- .useMart(
    biomart = biomart,
    dataset = dataset,
    host = host,
    verbose = FALSE,
    port = 443,
    ensemblRedirect = FALSE,
    http_config = http_config
  )

  return(ens)
}


## This function submits a small test query to identify a working Ensembl mirror.
## If no mirror argument is provided it will use "www" as its first choice.
## If the selected mirror returns a success (http 200) response it will be used
## Otherwise another mirror is selected at random and used instead.
## If all mirrors fail it will return an error
#' @importFrom httr2 req_body_form req_options req_timeout
#' @importFrom stringr str_match str_replace
.chooseEnsemblMirror <- function(mirror, http_config) {
  mirrors <- c("www", "asia", "useast")

  if (missing(http_config)) {
    http_config <- do.call(c, .getEnsemblSSL())
  }

  example_query <- '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "ensembl_gene_id" value = "ENSG00000000003"/>
		<Attribute name = "ensembl_gene_id" />
	</Dataset>
</Query>'

  ## create Ensembl URL and stop any redirection to a mirror
  host <- .constructEnsemblURL(mirror = mirror)
  host <- paste0(host, "/biomart/martservice?redirect=no")
  mirror <- str_match(host, pattern = "://([a-z]{3,6})\\.")[1, 2]

  req <- httr2::request(host) |>
    req_user_agent(
      .biomaRt_user_agent()
    ) |>
    req_body_form(query = example_query) |>
    req_timeout(10) |>
    req_options(!!!http_config)

  result <- tryCatch(httr2::req_perform(req), error = function(c) {
    "timeout"
  })

  tryAgain <- any(result == "timeout") || httr2::resp_status(result) == 500

  if (tryAgain) {
    ## try an alternative mirror if ensembl returns 500
    remaining_mirrors <- setdiff(mirrors, mirror)
    while ((length(remaining_mirrors) > 0) && (tryAgain)) {
      mirror <- sample(remaining_mirrors, size = 1)
      message("Ensembl site unresponsive, trying ", mirror, " mirror")
      host <- str_replace(
        host,
        pattern = "://([a-z]{3,6})\\.",
        replacement = paste0("://", mirror, ".")
      )

      req <- httr2::request(host) |>
        req_user_agent(
          .biomaRt_user_agent()
        ) |>
        req_body_form(query = example_query) |>
        req_timeout(10) |>
        req_options(!!!http_config)

      result <- tryCatch(httr2::req_perform(req), error = function(c) {
        "timeout"
      })
      tryAgain <- any(result == "timeout") || httr2::resp_status(result) == 500
      if (tryAgain) {
        remaining_mirrors <- setdiff(remaining_mirrors, mirror)
      }
    }
  }
  if (tryAgain) {
    stop("Unable to query any Ensembl site")
  }

  return(mirror)
}
