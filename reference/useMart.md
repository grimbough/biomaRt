# Connects to the selected BioMart database and dataset

A first step in using the biomaRt package is to select a BioMart
database and dataset to use. The useMart function enables one to connect
to a specified BioMart database and dataset within this database. To
know which BioMart databases are available see the
[`listMarts()`](https://huber-group-embl.github.io/biomaRt/reference/listMarts.md)
function. To know which datasets are available within a BioMart
database, first select the BioMart database using `useMart()` and then
use the
[`listDatasets()`](https://huber-group-embl.github.io/biomaRt/reference/listDatasets.md)
function on the selected BioMart, see
[`listDatasets()`](https://huber-group-embl.github.io/biomaRt/reference/listDatasets.md)
function.

## Usage

``` r
useMart(
  biomart,
  dataset,
  host = "https://www.ensembl.org",
  path = "/biomart/martservice",
  port,
  version,
  verbose = FALSE
)
```

## Arguments

- biomart:

  BioMart database name you want to connect to. Possible database names
  can be retrieved with the function
  [`listMarts()`](https://huber-group-embl.github.io/biomaRt/reference/listMarts.md)

- dataset:

  Dataset you want to use. To see the different datasets available
  within a biomaRt you can e.g. do: mart =
  [useMart()](https://huber-group-embl.github.io/biomaRt/reference/),
  followed by
  [listDatasets()](https://huber-group-embl.github.io/biomaRt/reference/mart).

- host:

  Host to connect to. Defaults to `www.ensembl.org`

- path:

  Path that should be pasted after to host to get access to the web
  service URL

- port:

  port to connect to, will be pasted between host and path

- version:

  Use version name instead of biomart name to specify which BioMart you
  want to use

- verbose:

  Give detailed output of what the method is doing while in use, for
  debugging

## Value

An object of class Mart, which can be used in functions such as
[`getBM()`](https://huber-group-embl.github.io/biomaRt/reference/getBM.md),
[`getLDS()`](https://huber-group-embl.github.io/biomaRt/reference/getLDS.md),
[`getGene()`](https://huber-group-embl.github.io/biomaRt/reference/getGene.md),
etc. to retrieve data from the selected BioMart database.

## Details

The previously available `archive` argument is defunct. A better
alternative is to specify the url of the archived BioMart you want to
access. For Ensembl you can view the list of archives using
[`listEnsemblArchives()`](https://huber-group-embl.github.io/biomaRt/reference/listEnsemblArchives.md).

## Author

Steffen Durinck, Mike L. Smith

## Examples

``` r
if (FALSE) { # interactive()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)
}
```
