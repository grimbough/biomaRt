# Gives a summary of the attribute pages

Attributes in BioMart databases are grouped together in attribute pages.
The `attributePages()` function gives a summary of the attribute
categories and groups present in the BioMart. These page names can be
used to display only a subset of the available attributes in the
[`listAttributes()`](https://huber-group-embl.github.io/biomaRt/reference/listAttributes.md)
function.

## Usage

``` r
attributePages(mart)
```

## Arguments

- mart:

  object of class Mart, created with the
  [`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md)
  function.

## Value

A character vector containing the names of the attribute pages present
in the `mart` object.

## Author

Steffen Durinck

## Examples

``` r
if (FALSE) { # interactive()
mart <- useMart(
  "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)
attributePages(mart)
}
```
