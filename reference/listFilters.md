# List or search the filters available in the selected dataset

Filters are what we use as inputs for a biomaRt query. For example, if
we want to retrieve all EntrezGene identifiers on chromosome X,
`chromosome` will be the filter, with corresponding value X.

## Usage

``` r
listFilters(mart, what = c("name", "description"))

searchFilters(mart, pattern = ".*")
```

## Arguments

- mart:

  object of class `Mart` created using the
  [`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md)
  function

- what:

  character vector indicating what information to display about the
  available filters. Valid values are `name`, `description`, `options`,
  `fullDescription`, `filters`, `type`, `operation`, `filters8`,
  `filters9`.

- pattern:

  Character vector defining the regular expression
  ([regex](https://rdrr.io/r/base/regex.html)) to be used for the
  search. If left blank the default is to use \`".\*"“ which will match
  everything.

## Author

Steffen Durinck, Mike Smith

## Examples

``` r
if (FALSE) { # interactive()
## list the available Ensembl marts and use Ensembl Genes
listEnsembl()
ensembl <- useEnsembl(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)

## list the available datasets in this Mart
listFilters(mart = ensembl)

## the list of filters is long and not easy to read
## we can search for a term of interest to reduce this e.g. 'gene'
searchFilters(mart = ensembl, pattern = "gene")

## search the available filters to find entries containing 'entrez' or 'hgnc'
searchFilters(mart = ensembl, 'entrez|hgnc')
}
```
