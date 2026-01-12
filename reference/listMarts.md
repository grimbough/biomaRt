# lists the available BioMart databases

This function returns a list of BioMart databases to which biomaRt can
connect. By default the Ensembl BioMart databases are displayed. To
establish a connection use the
[`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md)
function.

## Usage

``` r
listMarts(
  mart = NULL,
  host = "https://www.ensembl.org",
  path = "/biomart/martservice",
  port,
  includeHosts = FALSE,
  archive = FALSE,
  http_config = list(),
  verbose = FALSE
)
```

## Arguments

- mart:

  mart object created with the
  [`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md)
  function. This is optional, as you usually use `listMarts()` to see
  which marts there are to connect to.

- host:

  Host to connect to. Defaults to `www.ensembl.org`

- path:

  path to martservice that should be pasted behind the host to get to
  web service URL

- port:

  port to use in HTTP communication

- includeHosts:

  boolean to indicate if function should return host of the BioMart
  databases

- archive:

  Boolean to indicate if you want to access archived versions of BioMart
  database. Note that this argument is now defunct and setting this
  value to `TRUE` will produce an error. A better alternative is to
  specify the url of the archived BioMart you want to access. For
  Ensembl you can view the list of archives using
  [`listEnsemblArchives()`](https://huber-group-embl.github.io/biomaRt/reference/listEnsemblArchives.md)

- http_config:

  Some hosts require specific HTTP settings to be used when connecting.
  This argument takes the output of
  [`httr::config()`](https://httr.r-lib.org/reference/config.html) and
  will be used when connecting to `host`. Can be ignored if you
  experience no problems accessing `host`.

- verbose:

  Give detailed output of what the method is doing, for debugging
  purposes.

## Details

If you receive an error message saying 'Unexpected format to the list of
available marts', this is often because there is a problem with the
BioMart server you are trying to connect to, and something other than
the list of available marts is being returned - often some like a 'down
for maintenance' page. If you browse to the provided URL and find a page
that starts with '`<MartRegistry>`' this is the correct listing and you
should report the issue on the Bioconductor support site:
https://support.bioconductor.org

## Author

Steffen Durinck, Mike Smith

## Examples

``` r
if (FALSE) { # interactive()
listMarts()
}
```
