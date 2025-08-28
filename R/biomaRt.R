##########################
# biomaRt source code     #
##########################
#                        #
# Licence: Artistic       #
# Author: Steffen Durinck #
##########################

##############################################################
# martCheck                                                   #
#                                                            #
# This function checks if there is a valid Mart object,       #
# if a dataset is selected and                                #
# if the correct BioMart database has been selected (optional)#
##############################################################

martCheck <- function(mart, biomart = NULL) {
  if (missing(mart) || !inherits(mart, "Mart")) {
    stop(
      "You must provide a valid Mart object. To create a Mart object use the function: useMart.  Check ?useMart for more information."
    )
  }
  if (!is.null(biomart)) {
    martcheck <- martBM(mart)
    bmok <- FALSE
    for (k in seq_along(biomart)) {
      if (martcheck[1] == biomart[k]) {
        bmok <- TRUE
      }
    }
    if (!bmok) {
      stop(
        "This function only works when used with the ",
        biomart,
        " BioMart."
      )
    }
  }
  if (martDataset(mart) == "") {
    stop(
      "No dataset selected, please select a dataset first.  You can see the available datasets by using the listDatasets function see ?listDatasets for more information.  Then you should create the Mart object by using the useMart function.  See ?useMart for more information"
    )
  }
}


#' @importFrom httr2 req_options req_perform req_timeout resp_body_string
bmRequest <- function(request, http_config, verbose = FALSE) {
  if (verbose) {
    message("Attempting web service request:\n", request)
  }

  request <- httr2::request(request) |>
    req_timeout(getOption("timeout", default = 60)) |>
    req_options(!!!http_config)

  result <- req_perform(request)

  result2 <- resp_body_string(result)
  if (is.na(result2)) {
    result2 <- resp_body_string(result, encoding = "Latin1")
  }
  return(result2)
}

#######################################################
# listMarts:                                           #
# list all available BioMart databases by default      #
# listMarts will check the central service to see which#
# BioMart databases are present                        #
#######################################################

#' lists the available BioMart databases
#'
#' This function returns a list of BioMart databases to which biomaRt can
#' connect.  By default the Ensembl BioMart databases are displayed. To
#' establish a connection use the [useMart()] function.
#'
#' If you receive an error message saying 'Unexpected format to the list of
#' available marts', this is often because there is a problem with the BioMart
#' server you are trying to connect to, and something other than the list of
#' available marts is being returned - often some like a 'down for
#' maintainance' page.  If you browse to the provided URL and find a page that
#' starts with '`<MartRegistry>`' this is the correct listing and you
#' should report the issue on the Bioconductor support site:
#' https://support.bioconductor.org
#'
#' @param mart mart object created with the [useMart()] function.
#' This is optional, as you usually use [listMarts()] to see which
#' marts there are to connect to.
#' @param host Host to connect to. Defaults to `www.ensembl.org`
#' @param path path to martservice that should be pasted behind the host to get
#' to web service URL
#' @param port port to use in HTTP communication
#' @param includeHosts boolean to indicate if function should return host of
#' the BioMart databases
#' @param archive Boolean to indicate if you want to access archived versions
#' of BioMart database. Note that this argument is now defunct and setting this
#' value to `TRUE` will produce an error. A better alternative is to
#' specify the url of the archived BioMart you want to access.  For Ensembl you
#' can view the list of archives using [listEnsemblArchives()]
#' @param http_config Some hosts require specific HTTP settings to be used when
#' connecting. This argument takes the output of [httr::config()] and
#' will be used when connecting to `host`.  Can be ignored if you
#' experience no problems accessing `host`.
#' @param verbose Give detailed output of what the method is doing, for
#' debugging purposes.
#' @author Steffen Durinck, Mike Smith
#' @keywords methods
#' @examples
#'
#' if(interactive()){
#' listMarts()
#' }
#'
#' @export
listMarts <- function(
  mart = NULL,
  host = "https://www.ensembl.org",
  path = "/biomart/martservice",
  port,
  includeHosts = FALSE,
  archive = FALSE,
  http_config,
  verbose = FALSE
) {
  if (missing(port)) {
    port <- ifelse(startsWith(host, "https"), yes = 443, no = 80)
  }

  if (
    grepl(pattern = "^https://.*ensembl.org", x = host) &&
      missing(http_config)
  ) {
    http_config <- .getEnsemblSSL()
  }

  if (missing(http_config)) {
    http_config <- list()
  }

  .listMarts(
    mart = mart,
    host = host,
    path = path,
    port = port,
    includeHosts = includeHosts,
    archive = archive,
    verbose = verbose,
    http_config = http_config,
    ensemblRedirect = TRUE
  )
}

#' @importFrom methods is
.listMarts <- function(
  mart = NULL,
  host = "www.ensembl.org",
  path = "/biomart/martservice",
  port = 80,
  includeHosts = FALSE,
  archive = FALSE,
  verbose = FALSE,
  http_config,
  ensemblRedirect = NULL,
  warn = TRUE
) {
  request <- NULL
  if (is.null(mart)) {
    host <- .cleanHostURL(host, warn = warn)
    if (archive) {
      stop(
        "The archive = TRUE argument is now defunct.\n",
        "Use listEnsemblArchives() to find the URL to directly query an Ensembl archive."
      )
    }
    request <- paste0(
      host,
      ":",
      port,
      path,
      "?type=registry&requestid=biomaRt"
    )

    if (is(http_config, "list")) {
      http_config <- do.call(c, http_config)
    }
  } else if (is(mart, "Mart")) {
    request <- paste0(martHost(mart), "?type=registry&requestid=biomaRt")
    http_config <- martHTTPConfig(mart)
  } else {
    stop(
      mart,
      " object needs to be of class Mart created with the useMart function.\n",
      "If you don't have a Mart object yet, use listMarts() without arguments or only specify the host argument"
    )
  }

  if (!ensemblRedirect && grepl(x = request, pattern = "ensembl.org")) {
    request <- paste0(request, "&redirect=no")
  }

  registry <- bmRequest(
    request = request,
    http_config = http_config,
    verbose = verbose
  )

  ## check this looks like the MartRegistry XML, otherwise throw an error
  if (!grepl(x = registry, pattern = "^\n*<MartRegistry>")) {
    if (grepl(x = registry, pattern = "status.ensembl.org")) {
      stop(
        "Your query has been redirected to http://status.ensembl.org ",
        "indicating this Ensembl service is currently unavailable.",
        "\nLook at ?useEnsembl for details on how to try a mirror site.",
        call. = FALSE
      )
    } else {
      stop(
        "Unexpected format to the list of available marts.\n",
        "Please check the following URL manually, ",
        "and try ?listMarts for advice.\n",
        request,
        call. = FALSE
      )
    }
  }

  registry_xml2 <- xml2::read_xml(registry)
  registry_xml2 <- xml2::xml_children(registry_xml2)

  ## create a table with the registry information
  marts <- do.call("rbind", lapply(registry_xml2, FUN = xml2::xml_attrs))
  marts <- as.data.frame(marts[marts[, "visible"] == "1", , drop = FALSE])
  ## rename some columns
  names(marts)[names(marts) == "name"] <- "biomart"
  names(marts)[names(marts) == "displayName"] <- "version"
  names(marts)[names(marts) == "serverVirtualSchema"] <- "vschema"

  if (includeHosts) {
    return(as.list(marts))
  } else {
    return(marts[, c("biomart", "version")])
  }
}

#################################
# #                           # #
# # Generic BioMart functions # #
# #                           # #
#################################

#' Connects to the selected BioMart database and dataset
#'
#' A first step in using the biomaRt package is to select a BioMart database
#' and dataset to use.  The useMart function enables one to connect to a
#' specified BioMart database and dataset within this database.  To know which
#' BioMart databases are available see the [listMarts()] function.  To know which
#' datasets are available within a BioMart database, first select the BioMart
#' database using [useMart()] and then use the [listDatasets()] function on the
#' selected BioMart, see [listDatasets()] function.
#'
#'
#' @param biomart BioMart database name you want to connect to. Possible
#' database names can be retrieved with the functio [listMarts()]
#' @param dataset Dataset you want to use.  To see the different datasets
#' available within a biomaRt you can e.g. do: mart = [useMart()]('ensembl'),
#' followed by [listDatasets()](mart).
#' @param host Host to connect to. Defaults to `www.ensembl.org`
#' @param path Path that should be pasted after to host to get access to the
#' web service URL
#' @param port port to connect to, will be pasted between host and path
#' @param archive Boolean to indicate if you want to access archived versions
#' of BioMart databases.  Note that this argument is now deprecated and will be
#' removed in the future.  A better alternative is to leave archive = FALSE and
#' to specify the url of the archived BioMart you want to access.  For Ensembl
#' you can view the list of archives using [listEnsemblArchives()]
#' @param version Use version name instead of biomart name to specify which
#' BioMart you want to use
#' @param verbose Give detailed output of what the method is doing while in
#' use, for debugging
#' @author Steffen Durinck, Mike L. Smith
#' @keywords methods
#' @examples
#'
#' if(interactive()){
#'
#'     mart = useMart("ensembl")
#'     mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#' }
#'
#' @export
useMart <- function(
  biomart,
  dataset,
  host = "https://www.ensembl.org",
  path = "/biomart/martservice",
  port,
  archive = FALSE,
  version,
  verbose = FALSE
) {
  if (missing(port)) {
    port <- ifelse(startsWith(host, "https")[1], yes = 443, no = 80)
  }

  .useMart(
    biomart,
    dataset,
    host = host,
    path = path,
    port = port,
    archive = archive,
    version = version,
    verbose = verbose,
    http_config = list(),
    ensemblRedirect = TRUE
  )
}

.useMart <- function(
  biomart,
  dataset,
  host = "https://www.ensembl.org",
  path = "/biomart/martservice",
  port = 443,
  archive = FALSE,
  ensemblRedirect = NULL,
  version,
  http_config,
  verbose = FALSE
) {
  if (missing(biomart) && missing(version)) {
    stop(
      "No biomart databases specified. Specify a biomart database to use using the biomart or version argument"
    )
  }
  if (!missing(biomart) && !is.character(biomart)) {
    stop(
      "biomart argument is not a string. ",
      "The biomart argument should be a single character string"
    )
  }

  if (biomart == "ensembl" && grepl(x = host, pattern = "ensembl.org")) {
    biomart <- "ENSEMBL_MART_ENSEMBL"
  }

  reqHost <- host
  host <- .cleanHostURL(host)

  marts <- .listMarts(
    host = host,
    path = path,
    port = port,
    includeHosts = TRUE,
    http_config = http_config,
    archive = archive,
    ensemblRedirect = ensemblRedirect,
    warn = FALSE
  )
  mindex <- NA
  if (!missing(biomart)) {
    mindex <- match(biomart, marts$biomart)
  }
  if (!missing(version)) {
    mindex <- match(version, marts$version)
  }
  if (is.na(mindex) || archive) {
    mindex <- match(biomart, marts$database)
  }
  if (is.na(mindex)) {
    stop(
      "Incorrect BioMart name, use the listMarts function to see which BioMart databases are available"
    )
  }

  if (
    is.na(marts$path[mindex]) ||
      is.na(marts$vschema[mindex]) ||
      is.na(marts$host[mindex]) ||
      is.na(marts$port[mindex])
  ) {
    stop(
      "The selected biomart databases is not available due to error in the BioMart central registry, please report so the BioMart registry file can be fixed."
    )
  }

  if (marts$path[mindex] == "") {
    marts$path[mindex] <- "/biomart/martservice"
  } # temporary to catch bugs in registry

  if (!missing(version)) {
    biomart <- marts$biomart[mindex]
  }
  biomart <- sub(" ", "%20", biomart, fixed = TRUE, useBytes = TRUE)

  ## adding option to force use of specified host with ensembl
  redirect <- ifelse(
    !ensemblRedirect && grepl(x = host, pattern = "ensembl.org"),
    "?redirect=no",
    ""
  )

  if (missing(http_config)) {
    http_config <- list()
  }

  mart <- Mart(
    biomart = biomart,
    vschema = marts$vschema[mindex],
    host = paste0(host, ":", port, marts$path[mindex], redirect),
    http_config = http_config
  )

  if (any(grepl("archive", martHost(mart), fixed = TRUE))) {
    ## hack to work around redirection of most recent mirror URL
    archives <- .listEnsemblArchives(
      https = TRUE,
      http_config = http_config
    )
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

  BioMartVersion <- bmVersion(mart, verbose = verbose)

  if (verbose) {
    writeLines(paste(
      "BioMartServer running BioMart version:",
      BioMartVersion,
      sep = " "
    ))
    writeLines(paste("Mart virtual schema:", martVSchema(mart), sep = " "))
    if (!any(grepl(reqHost, martHost(mart), fixed = TRUE))) {
      writeLines(paste(
        "Requested host was redirected from ",
        reqHost,
        " to ",
        martHost(mart),
        sep = ""
      ))
    }
    writeLines(paste("Mart host:", martHost(mart), sep = " "))
  }
  if (!missing(dataset)) {
    mart <- useDataset(mart = mart, dataset = dataset, verbose = verbose)
  }
  return(mart)
}


#' List or search the datasets available in the selected BioMart database
#'
#' Lists or search the datasets available in the selected BioMart database
#'
#' @param mart object of class Mart created with the useMart function
#' @param verbose Give detailed output of what the method is doing, for
#' debugging purposes
#' @param pattern Character vector defining the regular expression
#' ([regex][base::regex]) to be used for the search.  If left blank the
#' default is to use ".*" which will match everything and return the same as
#' [listDatasets()].
#' @author Steffen Durinck, Mike Smith
#' @keywords methods
#' @examples
#'
#'
#' if(interactive()){
#'
#'     ## list the available Ensembl marts and use Ensembl Genes
#'     listEnsembl()
#'     ensembl <- useEnsembl(biomart = "ensembl")
#'
#'     ## list the available datasets in this Mart
#'     listDatasets(mart = ensembl)
#'
#'     ## the list of Ensembl datasets grows ever larger (101 as of Ensembl 93)
#'     ## we can search for a term of interest to reduce the length e.g. 'sapiens'
#'     searchDatasets(mart = ensembl, pattern = "sapiens")
#'
#'     ## search for any dataset containing the word Rat or rat
#'     searchDatasets(mart = ensembl, pattern = "(R|r)at")
#' }
#'
#' @export
listDatasets <- function(mart, verbose = FALSE) {
  .listDatasets(mart = mart, verbose = verbose, sort = TRUE)
}

#' @importFrom methods is
.listDatasets <- function(mart, verbose = FALSE, sort = FALSE) {
  if (missing(mart) || !is(mart, "Mart")) {
    stop("No Mart object given or object not of class 'Mart'")
  }

  ## we choose a separator based on whether 'redirect=no' is present
  ## should always be '?' now
  sep <- ifelse(grepl(x = martHost(mart), pattern = ".+\\?.+"), "&", "?")

  request <- paste0(
    martHost(mart),
    sep,
    "type=datasets&requestid=biomaRt&mart=",
    martBM(mart)
  )
  http_config <- martHTTPConfig(mart)

  bmResult <- bmRequest(
    request = request,
    http_config = http_config,
    verbose = verbose
  )
  con <- textConnection(bmResult)
  txt <- scan(
    con,
    sep = "\t",
    blank.lines.skip = TRUE,
    what = "character",
    quiet = TRUE,
    quote = "\""
  )
  close(con)

  ## select visible ("1") table sets
  i <- intersect(which(txt == "TableSet"), which(txt == "1") - 3L)

  res <- data.frame(
    dataset = I(txt[i + 1L]),
    description = I(txt[i + 2L]),
    version = I(txt[i + 4L])
  )

  ## sort alphabetically
  if (sort) {
    res <- res[order(res$dataset), ]
  }
  rownames(res) <- NULL

  return(res)
}

## Check version of BioMart service
#' @importFrom utils read.table
bmVersion <- function(mart, verbose = FALSE) {
  ## save some time and a HTTP request if this is Ensembl
  if (grepl(pattern = "ensembl.org", x = martHost(mart), fixed = TRUE)) {
    bmv <- "0.7"
  } else {
    ## we choose a separator based on whether 'redirect=no' is present
    sep <- ifelse(grepl(x = martHost(mart), pattern = ".+\\?.+"), "&", "?")

    request <- paste0(
      martHost(mart),
      sep,
      "type=version",
      "&requestid=biomaRt&mart=",
      martBM(mart)
    )
    http_config <- martHTTPConfig(mart)

    BioMartVersion <- bmRequest(
      request = request,
      http_config = http_config,
      verbose = verbose
    )
    bmv <- ""
    if (BioMartVersion == "\n" || BioMartVersion == "") {
      bmv <- NA
      if (verbose) {
        warning(paste(
          "BioMart version is not available from BioMart server:",
          request,
          sep = "\n"
        ))
      }
    } else {
      con <- textConnection(BioMartVersion)
      bmVersionParsed <- read.table(
        con,
        sep = "\t",
        header = FALSE,
        quote = "",
        comment.char = "",
        as.is = TRUE
      )
      close(con)
      if (verbose) {
        print(bmVersionParsed)
      }

      if (dim(bmVersionParsed)[2] >= 1) {
        bmv <- bmVersionParsed[1, 1]
      }
    }
  }
  return(bmv)
}


#' @importFrom utils read.table
.getAttrFilt <- function(mart, verbose, type) {
  ## we choose a separator based on whether 'redirect=no' is present
  sep <- ifelse(grepl(x = mart@host, pattern = ".+\\?.+"), "&", "?")

  request <- paste0(
    mart@host,
    sep,
    "type=",
    type,
    "&dataset=",
    martDataset(mart),
    "&requestid=biomaRt&mart=",
    martBM(mart),
    "&virtualSchema=",
    martVSchema(mart)
  )

  attrfilt <- bmRequest(
    request = request,
    http_config = martHTTPConfig(mart),
    verbose = verbose
  )
  attrfiltParsed <- read.table(
    text = attrfilt,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    as.is = TRUE
  )
  return(attrfiltParsed)
}

.getAttributes <- function(mart, verbose = FALSE) {
  attributes_table <- .getAttrFilt(
    mart = mart,
    verbose = verbose,
    type = "attributes"
  )

  if (ncol(attributes_table) < 4) {
    stop(
      "biomaRt error: looks like we're connecting to incompatible version of BioMart."
    )
  }

  colnames(attributes_table) <- c(
    "name",
    "description",
    "fullDescription",
    "page"
  )
  return(attributes_table)
}

.getFilters <- function(mart, verbose = FALSE) {
  filters_table <- .getAttrFilt(
    mart = mart,
    verbose = verbose,
    type = "filters"
  )

  if (ncol(filters_table) < 7) {
    stop(
      "biomaRt error: looks like we're connecting to incompatible version of BioMart."
    )
  }

  colnames(filters_table) <- c(
    "name",
    "description",
    "options",
    "fullDescription",
    "filters",
    "type",
    "operation"
  )
  return(filters_table)
}

## Utility function to check dataset specification
## Returns dataset name as a character assuming all checks
## have been passed.
checkDataset <- function(dataset, mart) {
  validDatasets <- .listDatasets(mart, sort = FALSE)
  ## subsetting data.frames can produce some weird classes
  ## which aren't character(), so we coerce it here
  dataset <- as.character(dataset)

  if (length(dataset) > 1) {
    stop("Please only specify a single dataset name")
  }

  if (is.na(match(dataset, validDatasets$dataset))) {
    stop(
      "The given dataset: ",
      dataset,
      ", is not valid.  Correct dataset names can be obtained with the listDatasets() function."
    )
  }

  return(dataset)
}

## Select a BioMart dataset

#' Select a dataset to use and updates Mart object
#'
#' This function selects a dataset and updates the Mart object
#'
#'
#' @param dataset Dataset you want to use.  List of possible datasets can be
#' retrieved using the function [listDatasets()]
#' @param mart Mart object created with the [useMart()] function
#' @param verbose Give detailed output of what the method is doing, for
#' debugging
#' @author Steffen Durinck
#' @keywords methods
#' @examples
#'
#' if(interactive()){
#' mart=useMart("ensembl")
#' mart=useDataset("hsapiens_gene_ensembl", mart = mart)
#' }
#'
#' @export
useDataset <- function(dataset, mart, verbose = FALSE) {
  if (missing(mart) || !inherits(mart, "Mart")) {
    stop(
      "No valid Mart object given, specify a Mart object with the attribute mart"
    )
  }

  if (missing(dataset)) {
    stop(
      "No dataset given.  Please use the dataset argument to specify which dataset you want to use. Correct dataset names can be obtained with the listDatasets() function."
    )
  }

  dataset <- checkDataset(dataset = dataset, mart = mart)
  martDataset(mart) <- dataset

  if (verbose) {
    message("Checking attributes ...", appendLF = FALSE)
  }
  martAttributes(mart) <- .getAttributes(mart, verbose = verbose)
  if (verbose) {
    message(" ok")
    message("Checking filters ...", appendLF = FALSE)
  }
  martFilters(mart) <- .getFilters(mart, verbose = verbose)
  if (verbose) {
    message(" ok")
  }
  return(mart)
}

## listAttributes

#' lists the attributes available in the selected dataset
#'
#' Attributes are the outputs of a biomaRt query, they are the information we
#' want to retrieve.  For example if we want to retrieve all EntrezGene
#' identifiers of genes located on chromosome X, `entrezgene_id` will be
#' the attribute we use in the query.  The `listAttributes` function lists
#' the available attributes in the selected dataset.
#'
#'
#' @param mart object of class Mart created using the [useMart()] function
#' @param page Show only the attributes that belong to the specified attribute
#' page.
#' @param what vector of types of information about the attributes that need to
#' be displayed.  Can have values like name, description, fullDescription, page
#' @param pattern Character vector defining the regular expression
#' ([regex][base::regex]) to be used for the search.  If left blank the
#' default is to use ".*" which will match everything.
#' @author Steffen Durinck, Mike Smith
#' @keywords methods
#' @examples
#'
#'
#' if(interactive()){
#'
#'     ## list the available Ensembl marts and use Ensembl Genes
#'     listEnsembl()
#'     ensembl <- useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')
#'
#'     ## list the available datasets in this Mart
#'     listAttributes(mart = ensembl)
#'
#'     ## the list of attributes is very long and gets truncated by R
#'     ## we can search for a term of interest to filter this e.g. 'start'
#'     searchAttributes(mart = ensembl, pattern = "start")
#'
#'     ## filter the attributes to give only entries containing 'entrez' or 'hgnc'
#'     searchAttributes(mart = ensembl, 'entrez|hgnc')
#' }
#'
#' @export
listAttributes <- function(
  mart,
  page,
  what = c("name", "description", "page")
) {
  martCheck(mart)
  if (!missing(page) && !page %in% attributePages(mart)) {
    stop(
      "The chosen page: ",
      page,
      " is not valid, please use the correct page name using the attributePages function"
    )
  }
  attrib <- NULL
  if (!missing(page)) {
    sel <- which(martAttributes(mart)[, "page"] == page)
    attrib <- martAttributes(mart)[sel, what]
  } else {
    attrib <- martAttributes(mart)[, what]
  }
  return(attrib)
}

## attributePages

#' Gives a summary of the attribute pages
#'
#' Attributes in BioMart databases are grouped together in attribute pages.
#' The [attributePages()] function gives a summary of the attribute categories and
#' groups present in the BioMart.  These page names can be used to display only
#' a subset of the available attributes in the [listAttributes()] function.
#'
#'
#' @param mart object of class Mart, created with the [useMart()] function.
#' @author Steffen Durinck
#' @keywords methods
#' @examples
#'
#'
#' if(interactive()){
#' mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
#' attributePages(mart)
#' }
#'
#' @export
attributePages <- function(mart) {
  martCheck(mart)
  pages <- unique(martAttributes(mart)[, "page"])
  return(pages)
}

## listFilters

#' List or search the filters available in the selected dataset
#'
#' Filters are what we use as inputs for a biomaRt query.  For example, if we
#' want to retrieve all EntrezGene identifiers on chromosome X,
#' `chromosome` will be the filter, with corresponding value X.
#'
#'
#' @param mart object of class `Mart` created using the
#' [useMart()] function
#' @param what character vector indicating what information to display about
#' the available filters.  Valid values are `name`, `description`,
#' `options`, `fullDescription`, `filters`, `type`,
#' `operation`, `filters8`, `filters9`.
#' @param pattern Character vector defining the regular expression
#' ([regex][base::regex]) to be used for the search.  If left blank the
#' default is to use `".*"`` which will match everything.
#' @author Steffen Durinck, Mike Smith
#' @keywords methods
#' @examples
#'
#'
#' if(interactive()){
#'
#'     ## list the available Ensembl marts and use Ensembl Genes
#'     listEnsembl()
#'     ensembl <- useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')
#'
#'     ## list the available datasets in this Mart
#'     listFilters(mart = ensembl)
#'
#'     ## the list of filters is long and not easy to read
#'     ## we can search for a term of interest to reduce this e.g. 'gene'
#'     searchFilters(mart = ensembl, pattern = "gene")
#'
#'     ## search the available filters to find entries containing 'entrez' or 'hgnc'
#'     searchFilters(mart = ensembl, 'entrez|hgnc')
#' }
#'
#' @export
listFilters <- function(mart, what = c("name", "description")) {
  martCheck(mart)
  filters <- martFilters(mart)
  badwhat <- !(what %in% colnames(filters))
  if (any(badwhat)) {
    stop(sprintf(
      "The function argument 'what' contains %s: %s\nValid are: %s\n",
      if (sum(badwhat) > 1) "invalid values" else "an invalid value",
      paste(what[badwhat], collapse = ", "),
      paste(colnames(filters), collapse = ", ")
    ))
  }
  return(filters[, what])
}

#' @export
filterOptions <- function(filter, mart) {
  .Defunct(
    new = "listFilterOptions",
    msg = c(
      "filterOptions() has been made defunct and will be removed from biomaRt.",
      "\nPlease use listFilterOptions() instead."
    )
  )
}

## filterType

#' Displays the filter type
#'
#' Displays the type of the filer given a filter name.
#'
#'
#' @param filter A valid filter name. Valid filters are given by the
#' [listFilters()] function
#' @param mart object of class Mart, created using the [useMart()] function
#' @author Steffen Durinck
#' @keywords methods
#' @examples
#'
#'
#' if(interactive()){
#' mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
#' filterType("chromosome_name", mart)
#' }
#'
#' @export
filterType <- function(filter, mart) {
  if (missing(filter)) {
    stop(
      "No filter given. Please specify the filter for which you want to retrieve the filter type"
    )
  }
  if (!is.character(filter)) {
    stop("Filter argument should be of class character")
  }
  martCheck(mart)
  type <- "unknown"
  sel <- which(listFilters(mart, what = "name") == filter)
  if (is.null(sel)) {
    stop(paste("Invalid filter", filter, sep = ": "))
  }
  type <- listFilters(mart, what = "type")[sel]
  return(type)
}

##########################################
# getBM: generic BioMart query function   #
##########################################

#' Retrieves information from the BioMart database
#'
#' This function is the main biomaRt query function.  Given a set of filters
#' and corresponding values, it retrieves the user specified attributes from
#' the BioMart database one is connected to.
#'
#'
#' @param attributes Attributes you want to retrieve.  A possible list of
#' attributes can be retrieved using the function [listAttributes()].
#' @param filters Filters (one or more) that should be used in the query.  A
#' possible list of filters can be retrieved using the function [listFilters()].
#' @param values Values of the filter, e.g. vector of affy IDs.  If multiple
#' filters are specified then the argument should be a list of vectors of which
#' the position of each vector corresponds to the position of the filters in
#' the filters argument.
#' @param mart object of class Mart, created with the [useMart()] function.
#' @param checkFilters Sometimes attributes where a value needs to be
#' specified, for example upstream_flank with value 20 for obtaining upstream
#' sequence flank regions of length 20bp, are treated as filters in BioMarts.
#' To enable such a query to work, one must specify the attribute as a filter
#' and set `checkFilters = FALSE` for the query to work.
#' @param verbose When using biomaRt in webservice mode and setting verbose to
#' TRUE, the XML query to the webservice will be printed.
#' @param uniqueRows If the result of a query contains multiple identical rows,
#' setting this argument to `TRUE` (default) will result in deleting the
#' duplicated rows in the query result at the server side.
#' @param bmHeader Boolean to indicate if the result retrieved from the BioMart
#' server should include the data headers or not, defaults to `FALSE`.  This
#' should only be switched on if the default behavior results in errors,
#' setting to on might still be able to retrieve your data in that case
#' @param quote Sometimes parsing of the results fails due to errors in the
#' Ensembl data fields such as containing a quote, in such cases you can try to
#' change the value of quote to try to still parse the results.
#' @param useCache Boolean indicating whether the results cache should be used.
#' Setting to `FALSE` will disable reading and writing of the cache.  This
#' argument is likely to disappear after the cache functionality has been
#' tested more thoroughly.
#' @return A `data.frame`. There is no implicit mapping between its rows
#' and the function arguments (e.g. `filters`, `values`), therefore
#' make sure to have the relevant identifier(s) returned by specifying them in
#' `attributes`. See Examples.
#' @author Steffen Durinck
#' @keywords methods
#' @examples
#'
#' if(interactive()){
#'   mart <- useEnsembl(biomart = "ensembl",
#'                      dataset = "hsapiens_gene_ensembl")
#'
#'   getBM(attributes = c("affy_hg_u95av2", "hgnc_symbol", "chromosome_name", "band"),
#'         filters    = "affy_hg_u95av2",
#'         values     = c("1939_at","1503_at","1454_at"),
#'         mart       = mart)
#'   }
#'
#' @importFrom progress progress_bar
#' @export
getBM <- function(
  attributes,
  filters = "",
  values = "",
  mart,
  checkFilters = TRUE,
  verbose = FALSE,
  uniqueRows = TRUE,
  bmHeader = FALSE,
  quote = "\"",
  useCache = TRUE
) {
  ## check the arguments are all valid
  martCheck(mart)
  if (missing(attributes)) {
    stop("Argument 'attributes' must be specified.")
  }

  if (is.list(filters) && !missing(values)) {
    warning(
      "Argument 'values' should not be used when argument 'filters' is a list and will be ignored."
    )
  }
  if (is.list(filters) && is.null(names(filters))) {
    stop("Argument 'filters' must be a named list when sent as a list.")
  }
  if (!is.list(filters) && all(nzchar(filters)) && missing(values)) {
    stop("Argument 'values' must be specified.")
  }
  if (length(filters) > 0 && length(values) == 0) {
    stop("Values argument contains no data.")
  }
  if (is.list(filters)) {
    values <- filters
    filters <- names(filters)
  }
  if (!is.logical(uniqueRows)) {
    stop(
      "Argument 'uniqueRows' must be a logical value, so either TRUE or FALSE"
    )
  }

  ## determine if we should use the results cache
  if (useCache) {
    cache <- .biomartCacheLocation()
    bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  }
  hash <- .createHash(mart, attributes, filters, values, uniqueRows, bmHeader)
  if (useCache && .checkValidCache(bfc, hash)) {
    if (verbose) {
      message("Cache found")
    }
    result <- .readFromCache(bfc, hash)
    return(result)
  }
  ## force the query to return the 'descriptive text' header names with the result
  ## we use these later to match and order attribute/column names
  xmlQuery <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "',
    martVSchema(mart),
    '" uniqueRows="',
    as.numeric(uniqueRows),
    '" count="0" datasetConfigVersion="0.6" header="1"',
    ' formatter="TSV" requestid="biomaRt"><Dataset name="',
    martDataset(mart),
    '">'
  )

  # checking the Attributes
  invalid <- !(attributes %in% listAttributes(mart, what = "name"))
  if (any(invalid)) {
    stop(
      "Invalid attribute(s):",
      paste(attributes[invalid], collapse = ", "),
      "\nPlease use the function 'listAttributes' to get valid attribute names"
    )
  }

  # attribute are ok lets add them to the query
  attributeXML <- paste0(
    '<Attribute name = "',
    attributes,
    '"/>',
    collapse = ""
  )

  # checking the filters
  if (filters[1] != "" && checkFilters) {
    invalid <- !(filters %in% listFilters(mart, what = "name"))
    if (any(invalid)) {
      stop(
        "Invalid filters(s): ",
        paste(filters[invalid], collapse = ", "),
        "\nPlease use the function 'listFilters' to get valid filter names"
      )
    }
  }

  ## filterXML is a list containing filters with reduced numbers of values
  ## to meet the 500 value limit in BioMart queries
  filterXmlList <- .generateFilterXML(filters, values, mart)

  resultList <- list()
  if (length(filterXmlList) > 1) {
    pb <- progress_bar$new(
      total = length(filterXmlList),
      width = getOption("width") - 10,
      format = "Batch submitting query [:bar] :percent eta: :eta"
    )
    pb$tick(0)
    on.exit(pb$terminate())
  }

  ## we submit a query for each chunk of the filter list
  for (i in seq_along(filterXmlList)) {
    if (i > 1) {
      pb$tick()
    }

    filterXML <- filterXmlList[[i]]
    fullXmlQuery <- paste(
      xmlQuery,
      attributeXML,
      filterXML,
      "</Dataset></Query>",
      sep = ""
    )

    if (verbose) {
      message(fullXmlQuery)
    }

    ## we choose a separator based on whether '?redirect=no' is present
    sep <- ifelse(
      grepl(x = martHost(mart), pattern = ".+\\?.+"),
      "&",
      "?"
    )

    ## create a unique name for this chunk & see if it has been run before
    chunk_hash <- tools::md5sum(
      bytes = charToRaw(paste(martHost(mart), fullXmlQuery))
    )
    tf <- file.path(
      tempdir(),
      paste0("biomaRt_tmp_", chunk_hash, ".rds")
    )
    if (!file.exists(tf)) {
      postRes <- .submitQueryXML(
        host = paste0(martHost(mart), sep),
        query = fullXmlQuery,
        http_config = martHTTPConfig(mart)
      )
      result <- .processResults(
        postRes,
        mart = mart,
        hostURLsep = sep,
        fullXmlQuery = fullXmlQuery,
        quote = quote,
        numAttributes = length(attributes)
      )
      saveRDS(result, file = tf)
    } else {
      result <- readRDS(tf)
    }
    resultList[[i]] <- .setResultColNames(
      result,
      mart = mart,
      attributes = attributes,
      bmHeader = bmHeader
    )
  }
  ## collate results
  result <- do.call("rbind", resultList)

  if (useCache) {
    .addToCache(bfc = bfc, result = result, hash = hash)
  }

  ## remove any temp chunk files
  file.remove(list.files(
    tempdir(),
    pattern = "^biomaRt.*rds$",
    full.names = TRUE
  ))
  return(result)
}

###################################
# getLDS: Multiple dataset linking #
###################################

#' Retrieves information from two linked datasets
#'
#' This function is the main biomaRt query function that links 2 datasets and
#' retrieves information from these linked BioMart datasets.  In Ensembl this
#' translates to homology mapping.
#'
#'
#' @param attributes Attributes you want to retrieve of primary dataset.  A
#' possible list of attributes can be retrieved using the function
#' [listAttributes()].
#' @param filters Filters that should be used in the query. These filters will
#' be applied to primary dataset.  A possible list of filters can be retrieved
#' using the function [listFilters()].
#' @param values Values of the filter, e.g. list of affy IDs
#' @param mart object of class Mart created with the [useMart()] function.
#' @param attributesL Attributes of linked dataset that needs to be retrieved
#' @param filtersL Filters to be applied to the linked dataset
#' @param valuesL Values for the linked dataset filters
#' @param martL Mart object representing linked dataset
#' @param verbose When using \pkg{biomaRt} in webservice mode and setting
#' verbose to `TRUE`, the XML query to the webservice will be printed.
#' Alternatively in MySQL mode the MySQL query will be printed.
#' @param uniqueRows Logical to indicate if the BioMart web service should
#' return unique rows only or not.  Has the value of either `TRUE` or `FALSE`
#' @param bmHeader Boolean to indicate if the result retrieved from the BioMart
#' server should include the data headers or not, defaults to `TRUE`.  This
#' should only be switched off if the default behavior results in errors,
#' setting to off might still be able to retrieve your data in that case
#' @author Steffen Durinck
#' @keywords methods
#' @examples
#'
#' if(interactive()){
#' human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#' mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#' getLDS(attributes = c("hgnc_symbol","chromosome_name", "start_position"),
#'     filters = "hgnc_symbol", values = "TP53", mart = human,
#'     attributesL = c("chromosome_name","start_position"), martL = mouse)
#' }
#'
#' @export
#' @importFrom methods is
#' @importFrom utils head read.table
getLDS <- function(
  attributes,
  filters = "",
  values = "",
  mart,
  attributesL,
  filtersL = "",
  valuesL = "",
  martL,
  verbose = FALSE,
  uniqueRows = TRUE,
  bmHeader = TRUE
) {
  martCheck(mart)
  martCheck(martL)

  if (martHost(mart) != martHost(martL)) {
    stop("Both datasets must be located on the same host.")
  }

  if (martBM(mart) != martBM(martL)) {
    stop(
      "Both datasets must be located in the same Mart.\n",
      "You are trying to combine datasets in ",
      martBM(mart),
      " and ",
      martBM(martL)
    )
  }

  invalid <- !(attributes %in% listAttributes(mart, what = "name"))
  if (any(invalid)) {
    stop(
      "Invalid attribute(s): ",
      paste(attributes[invalid], collapse = ", "),
      "\nPlease use the function 'listAttributes' to get valid attribute names"
    )
  }

  invalid <- !(attributesL %in% listAttributes(martL, what = "name"))
  if (any(invalid)) {
    stop(
      "Invalid attribute(s): ",
      paste(attributesL[invalid], collapse = ", "),
      "\nPlease use the function 'listAttributes' to get valid attribute names"
    )
  }

  if (filters[1] != "") {
    invalid <- !(filters %in% listFilters(mart, what = "name"))
    if (any(invalid)) {
      stop(
        "Invalid filters(s): ",
        paste(filters[invalid], collapse = ", "),
        "\nPlease use the function 'listFilters' to get valid filter names"
      )
    }
  }
  if (filtersL[1] != "") {
    invalid <- !(filtersL %in% listFilters(martL, what = "name"))
    if (any(invalid)) {
      stop(
        "Invalid filters(s): ",
        paste(filtersL[invalid], collapse = ", "),
        "\nPlease use the function 'listFilters' to get valid filter names"
      )
    }
  }

  xmlQuery <- sprintf(
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "%s" uniqueRows = "%s" count = "0" datasetConfigVersion = "0.6" header="%s" formatter = "TSV" requestid="biomaRt"> <Dataset name = "%s">',
    martVSchema(mart),
    as.numeric(uniqueRows),
    as.numeric(bmHeader),
    martDataset(mart)
  )

  attributeXML <- paste0(
    '<Attribute name = "',
    attributes,
    '"/>',
    collapse = ""
  )

  ## ignore the chunk size here
  filterXML <- .generateFilterXML(
    filters = filters,
    values = values,
    mart = mart,
    maxChunkSize = Inf
  )

  xmlQuery <- paste0(xmlQuery, attributeXML, filterXML, "</Dataset>")

  xmlQuery <- paste0(xmlQuery, '<Dataset name = "', martDataset(martL), '" >')
  linkedAttributeXML <- paste0(
    '<Attribute name = "',
    attributesL,
    '"/>',
    collapse = ""
  )
  linkedFilterXML <- .generateFilterXML(
    filters = filtersL,
    values = valuesL,
    mart = mart,
    maxChunkSize = Inf
  )

  xmlQuery <- paste0(
    xmlQuery,
    linkedAttributeXML,
    linkedFilterXML,
    "</Dataset></Query>"
  )

  if (verbose) {
    message(xmlQuery)
  }

  ## we choose a separator based on whether '?redirect=no' is present
  sep <- ifelse(grepl(x = martHost(mart), pattern = ".+\\?.+"), "&", "?")
  ## POST query
  postRes <- .submitQueryXML(
    host = paste0(martHost(mart), sep),
    query = xmlQuery,
    http_config = martHTTPConfig(mart)
  )

  if (any(startsWith(postRes, "^Query ERROR"))) {
    stop(postRes)
  }

  if (postRes != "") {
    con <- textConnection(postRes)
    result <- read.table(
      con,
      sep = "\t",
      header = bmHeader,
      quote = "\"",
      comment.char = "",
      as.is = TRUE,
      check.names = TRUE
    )
    close(con)

    if (nrow(result) > 0 && all(is.na(result[, ncol(result)]))) {
      result <- result[, -ncol(result), drop = FALSE]
    }

    res_attributes <- c(attributes, attributesL)
    if (
      !(is(result, "data.frame") &&
        (ncol(result) == length(res_attributes)))
    ) {
      print(head(result))
      stop(
        "The query to the BioMart webservice returned an invalid result: ",
        "the number of columns in the result table does not equal the number of attributes in the query. \n",
        "Please report this on the support site at http://support.bioconductor.org"
      )
    }
    if (!bmHeader) {
      # assumes order of results same as order of attibutes in input
      colnames(result) <- res_attributes
    }
  } else {
    warning("getLDS returns NULL.")
    result <- NULL
  }
  return(result)
}

####################
# export FASTA      #
####################

#' Exports getSequence results to FASTA format
#'
#' Exports getSequence results to FASTA format
#'
#'
#' @param sequences A data.frame that was the output of the [getSequence()]
#' function
#' @param file File to which you want to write the data
#' @author Steffen Durinck
#' @keywords methods
#' @examples
#'
#'
#' if(interactive()){
#'     mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#'
#'     #seq<-getSequence(chromosome=c(2,2),start=c(100000,30000),end=c(100300,30500),mart=mart)
#'     #exportFASTA(seq,file="test.fasta")
#'
#' }
#'
#' @export
exportFASTA <- function(sequences, file) {
  if (missing(sequences) || !is.data.frame(sequences)) {
    stop(
      "No data.frame given to write FASTA.  The data.frame should be the output of the getSequence function."
    )
  }
  if (missing(file)) {
    stop("Please provide filename to write to")
  }
  if (length(sequences[1, ]) == 2) {
    for (i in seq(along = sequences[, 2])) {
      cat(
        paste(">", sequences[i, 2], "\n", sep = ""),
        file = file,
        append = TRUE
      )
      cat(as.character(sequences[i, 1]), file = file, append = TRUE)
      cat("\n\n", file = file, append = TRUE)
    }
  } else {
    for (i in seq(along = sequences[, 2])) {
      cat(
        paste(
          ">chromosome_",
          sequences[i, 1],
          "_start_",
          sequences[i, 2],
          "_end_",
          sequences[i, 3],
          "\n",
          sep = ""
        ),
        file = file,
        append = TRUE
      )
      cat(as.character(sequences[i, 4]), file = file, append = TRUE)
      cat("\n\n", file = file, append = TRUE)
    }
  }
}

###################
# Nature Protocol
###################

#' Display the analysis code from the 2009 Nature protocols paper
#'
#' This function opens an editor displaying the analysis code of the Nature
#' Protocols 2009 paper
#'
#' The [edit()] function uses `getOption("editor")` to select
#' the editor. Use, for instance, `options(editor="emacs")` to set another
#' editor.
#'
#' @author Steffen Durinck, Wolfgang Huber
#' @seealso [edit()]
#' @keywords methods
#' @examples
#'
#' if(interactive()){
#' NP2009code()
#' }
#'
#' @export
#' @importFrom utils edit
NP2009code <- function() {
  edit(file = system.file("scripts", "Integration-NP.R", package = "biomaRt"))
}
