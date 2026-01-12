# Displays the filter type

Displays the type of the filer given a filter name.

## Usage

``` r
filterType(filter, mart)
```

## Arguments

- filter:

  A valid filter name. Valid filters are given by the
  [`listFilters()`](https://huber-group-embl.github.io/biomaRt/reference/listFilters.md)
  function

- mart:

  object of class Mart, created using the
  [`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md)
  function

## Author

Steffen Durinck

## Examples

``` r
if (FALSE) { # interactive()
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
filterType("chromosome_name", mart)
}
```
