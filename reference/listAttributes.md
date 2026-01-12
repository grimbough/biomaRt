# lists the attributes available in the selected dataset

Attributes are the outputs of a biomaRt query, they are the information
we want to retrieve. For example if we want to retrieve all EntrezGene
identifiers of genes located on chromosome X, `entrezgene_id` will be
the attribute we use in the query. The `listAttributes` function lists
the available attributes in the selected dataset.

## Usage

``` r
listAttributes(mart, page, what = c("name", "description", "page"))

searchAttributes(mart, pattern = ".*")
```

## Arguments

- mart:

  object of class Mart created using the
  [`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md)
  function

- page:

  Show only the attributes that belong to the specified attribute page.

- what:

  vector of types of information about the attributes that need to be
  displayed. Can have values like name, description, fullDescription,
  page

- pattern:

  Character vector defining the regular expression
  ([regex](https://rdrr.io/r/base/regex.html)) to be used for the
  search. If left blank the default is to use ".\*" which will match
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
  dataset = 'hsapiens_gene_ensembl'
)

## list the available datasets in this Mart
listAttributes(mart = ensembl)

## the list of attributes is very long and gets truncated by R
## we can search for a term of interest to filter this e.g. 'start'
searchAttributes(mart = ensembl, pattern = "start")

## filter the attributes to give only entries containing 'entrez' or 'hgnc'
searchAttributes(mart = ensembl, 'entrez|hgnc')
}
```
