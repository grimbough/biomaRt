## sometimes results can be returned by getBM() in a different order to we
## asked for them, which messes up the column names.  Here we try to match
## results to known attribute names and rename accordingly.
.setResultColNames <- function(result, mart, attributes, bmHeader = FALSE) {
  ## get all available attributes and
  ## filter only for the ones we've actually asked for
  att <- listAttributes(mart, what = c("name", "description"))
  att <- att[which(att[, "name"] %in% attributes), ]
  if (
    length(which(duplicated(att[, "description"]))) >
      length(which(duplicated(att)))
  ) {
    warning(
      "Cannot unambiguously match attribute names
                Ignoring bmHeader argument and using biomart
                description field"
    )
    return(result)
  }

  resultNames <- colnames(result)
  ## match the returned column names with the attribute names
  matches <- match(resultNames, att[, 2], NA)
  if (anyNA(matches)) {
    warning(
      "Problems assigning column names.",
      "Currently using the biomart description field.",
      "You may wish to set these manually."
    )
    return(result)
  }
  ## if we want to use the attribute names we specified, do this,
  ## otherwise we use the header returned with the query
  if (!bmHeader) {
    colnames(result) <- att[matches, 1]
  }
  ## now put things in the order we actually asked for the attributes in
  result <- result[, match(att[matches, 1], attributes), drop = FALSE]

  return(result)
}

## BioMart doesn't work well if the list of values provided to a filter is
## longer than 500 values.  It returns only a subset of the requested data
## and does so silently!  This function is designed to take a list of provided
## filters, and split any longer than 'maxChunkSize'.  It operates recursively
## incase there are multiple filters that need splitting, and should ensure
## all possible groupings of filters are retained.
.splitValues <- function(valuesList, maxChunkSize = 500) {
  vLength <- lengths(valuesList[[1]])

  if (all(vLength <= maxChunkSize)) {
    return(valuesList)
  }
  ## pick the next filter to split
  vIdx <- min(which(vLength > maxChunkSize))

  nchunks <- (vLength[vIdx] %/% maxChunkSize) + 1
  splitIdx <- rep(1:nchunks, each = ceiling(vLength[vIdx] / nchunks))[
    1:vLength[vIdx]
  ]

  ## a new list we will populate with the chunks
  tmpList <- list()
  for (i in 1:nchunks) {
    for (j in seq_along(valuesList)) {
      listIdx <- ((i - 1) * length(valuesList)) + j
      tmpList[[listIdx]] <- valuesList[[j]]
      tmpList[[listIdx]][[vIdx]] <- tmpList[[listIdx]][[vIdx]][which(
        splitIdx == i
      )]
    }
  }
  ## recursively call the function to process next filter
  valuesList <- .splitValues(tmpList, maxChunkSize = maxChunkSize)
  return(valuesList)
}

## Creating the filter XML for a single chunk of values.  Returns a character
## vector containing the XML lines for all specified filters & their
## attributes spliced together into a single string.
.createFilterXMLchunk <- function(filterChunk, mart) {
  individualFilters <- vapply(
    names(filterChunk),
    FUN = function(filter, values, mart) {
      ## if the filter exists and is boolean we do this
      if (
        filter %in%
          listFilters(mart, what = "name") &&
          grepl(
            "boolean",
            filterType(filter = filter, mart = mart),
            fixed = TRUE
          )
      ) {
        if (!is.logical(values[[filter]])) {
          stop(
            "'",
            filter,
            "' is a boolean filter and needs a ",
            "corresponding logical value of TRUE or FALSE to ",
            "indicate if the query should retrieve all data that ",
            "fulfill the boolean or alternatively that all data ",
            "that not fulfill the requirement should be retrieved.",
            call. = FALSE
          )
        }
        val <- as.numeric(!values[[filter]])
        val <- paste0("\" excluded = \"", val, "\" ")
      } else {
        ## otherwise the filter isn't boolean, or doesn't exist

        if (is.numeric(values[[filter]])) {
          values[[filter]] <- as.integer(values[[filter]])
        }
        val <- paste0(values[[filter]], collapse = ",")
        ## convert " to ' to avoid truncating the query string
        val <- gsub(x = val, pattern = "\"", replacement = "'", fixed = TRUE)
        val <- paste0('" value = "', val, '" ')
      }
      filterXML <- paste0("<Filter name = \"", filter, val, "/>")
      return(filterXML)
    },
    FUN.VALUE = character(1),
    filterChunk,
    mart,
    USE.NAMES = FALSE
  )

  filterXML <- paste0(individualFilters, collapse = "")
  return(filterXML)
}

.generateFilterXML <- function(
  filters = "",
  values,
  mart,
  maxChunkSize = 5000
) {
  ## return empty string if no filter specified & this isn't ensembl
  ## specifying no filter is generally bad, as it will get 'everything'
  ## and we might encounter the time out problem
  if (filters[1] == "") {
    return("")
  }
  ## if we have multiple filters, the values must be specified as a list.
  if (length(filters) > 1 && !is.list(values)) {
    stop(
      "If using multiple filters, the 'value' has to be a list.",
      "\nFor example, a valid list for 'value' could be: list(affyid=c('1939_at','1000_at'), chromosome= '16')",
      "\nHere we select on Affymetrix identifier and chromosome, only results that pass both filters will be returned"
    )
  }
  ## it's easy to not realise you're passing a data frame here, so check
  if (is.data.frame(values) && ncol(values) == 1) {
    values <- values[, 1]
  }

  if (!is.list(values)) {
    values <- list(values)
  }
  names(values) <- filters

  values <- .splitValues(list(values), maxChunkSize = maxChunkSize)

  filterXML_list <- lapply(values, .createFilterXMLchunk, mart)

  return(filterXML_list)
}

## it seems like pretty common practice for users to copy and paste the host
## name from a browser if they're not accessing Ensembl.  Typically this will
## include the "http://" and maybe a trailing "/" and this messes up our
## paste the complete URL strategy and produces something invalid.
## This function tidies that up to catch common variants.
.cleanHostURL <- function(host) {
  if (!grepl("^https?://", x = host)) {
    host <- paste0("http://", host)
  }

  parsed_url <- httr2::url_parse(host)

  ## just supplying 'ensembl.org' is no longer handled correctly
  ## stick 'www' in front if we see this
  if (parsed_url$hostname == "ensembl.org") {
    parsed_url$hostname <- "www.ensembl.org"
  }

  ## For HTTPS on Ensembl
  if (
    grepl("ensembl", parsed_url$hostname, fixed = TRUE) &&
      parsed_url$scheme != "https"
  ) {
    parsed_url$scheme <- "https"
  }

  host <- httr2::url_build(parsed_url)

  ## strip trailing slash
  host <- gsub(pattern = "/$", replacement = "", x = host)
  return(host)
}

.createErrorMessage <- function(error_code, host = "") {
  ## if we encounter internal server error, suggest using a mirror
  if (error_code == 500) {
    err_msg <- "biomaRt has encountered an unexpected server error."
  } else if (error_code == 509) {
    err_msg <- "biomaRt has exceeded the bandwidth allowance with this server."
  } else {
    err_msg <- paste0(
      "biomaRt has encountered an unknown server error. HTTP error code: ",
      error_code,
      "\nPlease report this on the Bioconductor support site at https://support.bioconductor.org/"
    )
  }

  if (grepl("ensembl", x = host, fixed = TRUE)) {
    err_msg <- c(
      err_msg,
      "\nConsider trying one of the Ensembl mirrors (for more details look at ?useEnsembl)"
    )
  }

  return(err_msg)
}

#' @importFrom httr2 req_body_form req_options req_timeout resp_body_string resp_status
.submitQueryXML <- function(host, query, http_config) {
  req <- httr2::request(host) |>
    req_body_form(query = query) |>
    req_timeout(max(getOption("timeout", default = 300), 300)) |>
    req_options(!!!http_config)

  res <- httr2::req_perform(req)

  if (httr2::resp_url(res) != host) {
    req2 <- req |>
      httr2::req_url(httr2::resp_url(res))
    res <- httr2::req_perform(req2)
  }

  if (httr2::resp_is_error(res)) {
    err_msg <- .createErrorMessage(error_code = resp_status(res), host = host)
    stop(err_msg, call. = FALSE)
  }

  ## content() prints a message about encoding not being supplied
  ## for ensembl.org - no default, so we suppress it
  return(resp_body_string(res))
}

## if parsing of TSV results fails, try this
.fetchHTMLresults <- function(host, query, http_config) {
  query <- gsub(x = query, pattern = "TSV", replacement = "HTML", fixed = TRUE)
  html_res <- .submitQueryXML(host, query, http_config)

  html <- xml2::read_html(html_res)
  table <- xml2::xml_find_first(html, ".//table")
  rows <- xml2::xml_find_all(table, ".//tr")
  cells <- lapply(rows, xml2::xml_find_all, ".//td|.//th")

  list_of_cells <- lapply(cells, FUN = xml2::xml_text)
  colnames <- list_of_cells[[1]]
  out <- as.data.frame(do.call(rbind, list_of_cells), row.names = NULL)
  out <- out[-1, ]
  colnames(out) <- colnames
  rownames(out) <- NULL
  return(out)
}

#' @param postRes Character vector of length 1 returned by server.  We expect
#' this to be a tab delimited string that comprises the whole table of results
#' including column headers.
#'
#' @noRd
#' @importFrom methods is
#' @importFrom utils read.table
.processResults <- function(
  postRes,
  mart,
  hostURLsep = "?",
  fullXmlQuery,
  quote = "\"",
  numAttributes
) {
  ## we expect only a character vector of length 1
  if (!is.character(postRes) || length(postRes) != 1L) {
    cli::cli_abort(
      c(
        "The query to the BioMart webservice returned an invalid result.",
        "i" = "biomaRt expected a character string of length 1.",
        "i" = "Please report this on the support site at
           {.url https://support.bioconductor.org}."
      )
    )
  }

  if (startsWith(postRes, "Query ERROR")) {
    stop(postRes)
  }

  ## convert the serialized table into a dataframe
  result <- tryCatch(
    read.table(
      text = postRes,
      sep = "\t",
      header = TRUE,
      quote = quote,
      comment.char = "",
      tryLogical = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) {
      ## if the error relates to number of element, try reading HTML version
      if (!grepl(x = e, pattern = "line [0-9]+ did not have [0-9]+ elements")) {
        stop(e)
      }
      .fetchHTMLresults(
        host = paste0(martHost(mart), hostURLsep),
        query = fullXmlQuery,
        http_config = martHTTPConfig(mart)
      )
    }
  )

  if (!is.data.frame(result) || ncol(result) != numAttributes) {
    cli::cli_abort(
      c(
        "The query to the BioMart webservice returned an invalid result.",
        "i" = "The number of columns in the result table does not equal the
              number of attributes in the query.",
        "i" = "Please report this on the support site at
              {.url https://support.bioconductor.org}."
      )
    )
  }

  return(result)
}

##############################################
## searching Attributes, Filters, and Datasets
##############################################

## given a data.frame, searches every column for
## the value in 'pattern'
## returns index of rows containing a match
.searchInternal <- function(pattern, data) {
  colIdx <- vapply(
    data,
    FUN = stringr::str_detect,
    FUN.VALUE = logical(length = nrow(data)),
    pattern = pattern
  )
  rowIdx <- apply(colIdx, 1, any)

  ## return either the matching rows, or NULL
  if (any(rowIdx)) {
    return(data[rowIdx, ])
  } else {
    cli::cli_inform("No matching datasets found.")
    return(NULL)
  }
}

#' @rdname listDatasets
#' @export
searchDatasets <- function(mart, pattern = ".*") {
  if (missing(mart)) {
    cli::cli_abort("Argument {.arg mart} must be specified.")
  }

  datasets <- listDatasets(mart)
  res <- .searchInternal(pattern = pattern, data = datasets)

  if (is.null(res)) {
    return(invisible(NULL))
  }
  res
}

#' @rdname listAttributes
#'
#' @export
searchAttributes <- function(mart, pattern = ".*") {
  if (missing(mart)) {
    cli::cli_abort("Argument {.arg mart} must be specified.")
  }

  attributes <- listAttributes(mart)
  res <- .searchInternal(pattern = pattern, data = attributes)

  if (is.null(res)) {
    return(invisible(NULL))
  }
  res
}

#' @rdname listFilters
#' @export
searchFilters <- function(mart, pattern = ".*") {
  if (missing(mart)) {
    cli::cli_abort("Argument {.arg mart} must be specified.")
  }

  filters <- listFilters(mart)
  res <- .searchInternal(pattern = pattern, data = filters)

  if (is.null(res)) {
    return(invisible(NULL))
  }
  res
}


## Some filters have a predefined list of options that can be selected.
## This function lets us search those values, given a specified filter.
#' @rdname listFilterOptions
#' @export
searchFilterOptions <- function(mart, filter, pattern = ".*") {
  if (missing(mart)) {
    cli::cli_abort("Argument {.arg mart} must be specified.")
  }
  if (missing(filter)) {
    cli::cli_abort("Argument {.arg filter} must be specified.")
  }

  ## first get all filters & their options, then reduce to what's requested
  filters <- listFilters(mart, what = c("name", "options"))
  filters <- filters[filters$name == filter, ]
  if (nrow(filters) == 0) {
    cli::cli_abort("Filter {.val {filter}} not found.")
  }
  options <- gsub(filters$options, pattern = "^\\[|\\]$", replacement = "")
  options <- strsplit(options, split = ",", fixed = TRUE)[[1]]

  res <- grep(x = options, pattern = pattern, ignore.case = TRUE, value = TRUE)

  if (length(res) == 0) {
    cli::cli_inform("No matching values found")
  } else {
    res
  }
}

#' @export
searchFilterValues <- function(mart, filter, pattern) {
  .Defunct(
    new = "listFilterOptions",
    msg = "This function has been renamed searchFilterOptions()"
  )
}


#' List or search the options available for a specified filter.
#'
#' Some filters have a predefined list of values that can be used to search
#' them.  These functions give access to this list of options for a named
#' filter, so you can check in the case where your biomaRt query is not finding
#' anything.
#'
#'
#' @param mart object of class `Mart` created using the
#' [useMart()], or [useEnsembl()] functions
#' @param filter The name of the filter whose options should be listed or
#' searched.  You can list available filters via [listFilters()]
#' @param pattern Character vector defining the regular expression
#' ([regex][base::regex]) to be used for the search.  If left blank the
#' default is to use ".*" which will match everything.
#' @author Mike Smith
#' @seealso [listFilters()]
#' @keywords methods
#'
#' @examplesIf interactive()
#' ## Use the Ensembl human genes dataset
#' ensembl <- useEnsembl(
#'   biomart = "ENSEMBL_MART_ENSEMBL",
#'   dataset = "hsapiens_gene_ensembl"
#' )
#'
#' ## we can search for the name of a filter we're interested in e.g. 'phenotype'
#' ## we need to use the name of the filter in the next function
#' searchFilters(ensembl, pattern = "phenotype")
#'
#' ## list all the options available to the 'phenotype_source' filter
#' listFilterOptions(mart = ensembl, filter = "phenotype_source")
#'
#' ## search the 'phenotype_description' filter for the term 'crohn'
#' searchFilterOptions(
#'   mart = ensembl,
#'   filter = "phenotype_description",
#'   pattern = "crohn"
#' )
#'
#' @export
listFilterOptions <- function(mart, filter) {
  searchFilterOptions(mart = mart, filter = filter)
}

#' @export
listFilterValues <- function(mart, filter) {
  .Defunct(
    new = "listFilterOptions",
    msg = "This function has been renamed listFilterOptions()"
  )
}
