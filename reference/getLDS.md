# Retrieves information from two linked datasets

This function is the main biomaRt query function that links 2 datasets
and retrieves information from these linked BioMart datasets. In Ensembl
this translates to homology mapping.

## Usage

``` r
getLDS(
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
)
```

## Arguments

- attributes:

  Attributes you want to retrieve of primary dataset. A possible list of
  attributes can be retrieved using the function
  [`listAttributes()`](https://huber-group-embl.github.io/biomaRt/reference/listAttributes.md).

- filters:

  Filters that should be used in the query. These filters will be
  applied to primary dataset. A possible list of filters can be
  retrieved using the function
  [`listFilters()`](https://huber-group-embl.github.io/biomaRt/reference/listFilters.md).

- values:

  Values of the filter, e.g. list of affy IDs

- mart:

  object of class Mart created with the
  [`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md)
  function.

- attributesL:

  Attributes of linked dataset that needs to be retrieved

- filtersL:

  Filters to be applied to the linked dataset

- valuesL:

  Values for the linked dataset filters

- martL:

  Mart object representing linked dataset

- verbose:

  When using biomaRt in webservice mode and setting verbose to `TRUE`,
  the XML query to the webservice will be printed. Alternatively in
  MySQL mode the MySQL query will be printed.

- uniqueRows:

  Logical to indicate if the BioMart web service should return unique
  rows only or not. Has the value of either `TRUE` or `FALSE`

- bmHeader:

  Boolean to indicate if the result retrieved from the BioMart server
  should include the data headers or not, defaults to `TRUE`. This
  should only be switched off if the default behavior results in errors,
  setting to off might still be able to retrieve your data in that case

## Author

Steffen Durinck

## Examples

``` r
if (FALSE) { # interactive()
human <- useMart(
  "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "https://dec2021.archive.ensembl.org"
)
mouse <- useMart(
  "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = "https://dec2021.archive.ensembl.org"
)
getLDS(
  attributes = c("hgnc_symbol","chromosome_name", "start_position"),
  filters = "hgnc_symbol",
  values = "TP53",
  mart = human,
  attributesL = c("chromosome_name","start_position"),
  martL = mouse
)
}
```
