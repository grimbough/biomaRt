# lists the available BioMart databases hosted by Ensembl

This function returns a list of BioMart databases hosted by Ensembl. To
establish a connection use the
[`useEnsembl()`](https://huber-group-embl.github.io/biomaRt/reference/useEnsembl.md)
function.

## Usage

``` r
listEnsembl(
  mart = NULL,
  version = NULL,
  GRCh = NULL,
  mirror = NULL,
  verbose = FALSE
)

listEnsemblGenomes(includeHosts = FALSE, host = NULL)
```

## Arguments

- mart:

  mart object created with the useEnsembl function. This is optional, as
  you usually use
  [`listMarts()`](https://huber-group-embl.github.io/biomaRt/reference/listMarts.md)
  to see which marts there are to connect to.

- version:

  Ensembl version to connect to when wanting to connect to an archived
  Ensembl version

- GRCh:

  GRCh version to connect to if not the current GRCh38, currently this
  can only be 37

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  'www', 'useast', 'asia'. If no mirror is specified the primary site at
  www.ensembl.org will be used.

- verbose:

  Give detailed output of what the method is doing, for debugging
  purposes

- includeHosts:

  If this option is set to `TRUE` a more detailed output is produced,
  including the URL used to access the corresponding mart.

- host:

  Host to connect to. Use this argument to specify and archive site for
  `listEnsemblGenomes()` to work with.

## Author

Steffen Durinck, Mike L. Smith

## Examples

``` r
if (FALSE) { # interactive()
listEnsembl()

## list the default Ensembl Genomes marts
listEnsemblGenomes()

## list only the marts available in the Ensmbl Plants 56 archive
listEnsemblGenomes(host = "https://eg56-plants.ensembl.org/")
}
```
