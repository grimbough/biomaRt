# List or search the datasets available in the selected BioMart database

Lists or search the datasets available in the selected BioMart database

## Usage

``` r
listDatasets(mart, verbose = FALSE)

searchDatasets(mart, pattern = ".*")
```

## Arguments

- mart:

  object of class Mart created with the useMart function

- verbose:

  Give detailed output of what the method is doing, for debugging
  purposes

- pattern:

  Character vector defining the regular expression
  ([regex](https://rdrr.io/r/base/regex.html)) to be used for the
  search. If left blank the default is to use ".\*" which will match
  everything and return the same as `listDatasets()`.

## Author

Steffen Durinck, Mike Smith

## Examples

``` r
if (FALSE) { # interactive()
## list the available Ensembl marts and use Ensembl Genes
listEnsembl()
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")

## list the available datasets in this Mart
listDatasets(mart = ensembl)

## the list of Ensembl datasets grows ever larger (101 as of Ensembl 93)
## we can search for a term of interest to reduce the length e.g. 'sapiens'
searchDatasets(mart = ensembl, pattern = "sapiens")

## search for any dataset containing the word Rat or rat
searchDatasets(mart = ensembl, pattern = "(R|r)at")
}
```
