# List or search the options available for a specified filter.

Some filters have a predefined list of values that can be used to search
them. These functions give access to this list of options for a named
filter, so you can check in the case where your biomaRt query is not
finding anything.

## Usage

``` r
searchFilterOptions(mart, filter, pattern = ".*")

listFilterOptions(mart, filter)
```

## Arguments

- mart:

  object of class `Mart` created using the
  [`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md),
  or
  [`useEnsembl()`](https://huber-group-embl.github.io/biomaRt/reference/useEnsembl.md)
  functions

- filter:

  The name of the filter whose options should be listed or searched. You
  can list available filters via
  [`listFilters()`](https://huber-group-embl.github.io/biomaRt/reference/listFilters.md)

- pattern:

  Character vector defining the regular expression
  ([regex](https://rdrr.io/r/base/regex.html)) to be used for the
  search. If left blank the default is to use ".\*" which will match
  everything.

## See also

[`listFilters()`](https://huber-group-embl.github.io/biomaRt/reference/listFilters.md)

## Author

Mike Smith

## Examples

``` r
if (FALSE) { # interactive()
## Use the Ensembl human genes dataset
ensembl <- useEnsembl(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)

## we can search for the name of a filter we're interested in e.g. 'phenotype'
## we need to use the name of the filter in the next function
searchFilters(ensembl, pattern = "phenotype")

## list all the options available to the 'phenotype_source' filter
listFilterOptions(mart = ensembl, filter = "phenotype_source")

## search the 'phenotype_description' filter for the term 'crohn'
searchFilterOptions(
  mart = ensembl,
  filter = "phenotype_description",
  pattern = "crohn"
)
}
```
