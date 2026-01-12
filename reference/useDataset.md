# Select a dataset to use and updates Mart object

This function selects a dataset and updates the Mart object

## Usage

``` r
useDataset(dataset, mart, verbose = FALSE)
```

## Arguments

- dataset:

  Dataset you want to use. List of possible datasets can be retrieved
  using the function
  [`listDatasets()`](https://huber-group-embl.github.io/biomaRt/reference/listDatasets.md)

- mart:

  Mart object created with the
  [`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md)
  function

- verbose:

  Give detailed output of what the method is doing, for debugging

## Author

Steffen Durinck

## Examples

``` r
if (FALSE) { # interactive()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
}
```
