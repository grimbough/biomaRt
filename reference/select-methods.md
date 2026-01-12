# Retrieve information from the BioMart databases

`select`, `columns` and `keys` are used together to extract data from a
`Mart` object. These functions work much the same as the classic biomaRt
functions such as
[`getBM()`](https://huber-group-embl.github.io/biomaRt/reference/getBM.md)
etc. and are provide here to make this easier for people who are
comfortable using these methods from other Annotation packages. Examples
of other objects in other packages where you can use these methods
include (but are not limited to): `ChipDb`, `OrgDb` `GODb`,
`InparanoidDb` and `ReactomeDb`.

## Usage

``` r
# S4 method for class 'Mart'
keys(x, keytype, ...)

# S4 method for class 'Mart'
keytypes(x)

# S4 method for class 'Mart'
columns(x)

# S4 method for class 'Mart'
select(x, keys, columns, keytype, ...)
```

## Arguments

- x:

  the `Mart` object. The dataset of the `Mart` object must already be
  specified for all of these methods.

- keytype:

  the keytype that matches the keys used. For the `select` methods, this
  is used to indicate the kind of ID being used with the keys argument.
  For the `keys` method this is used to indicate which kind of keys are
  desired from `keys`

- ...:

  other arguments. These include:

  pattern:

  :   the pattern to match (used by keys)

  column:

  :   the column to search on. This is used by keys and is for when the
      thing you want to pattern match is different from the keytype, or
      when you want to simply want to get keys that have a value for the
      thing specified by the column argument.

  fuzzy:

  :   TRUE or FALSE value. Use fuzzy matching? (this is used with
      pattern by the keys method)

- keys:

  the keys to select records for from the database. Keys for some
  keytypes can be extracted by using the `keys` method.

- columns:

  the columns or kinds of things that can be retrieved from the
  database. As with `keys`, all possible columns are returned by using
  the `columns` method.

## Value

`keys`,`columns` and `keytypes` each return a character vector or
possible values. `select` returns a data.frame.

## Details

`columns` shows which kinds of data can be returned from the `Mart`
object.

`keytypes` allows the user to discover which keytypes can be passed in
to `select` or `keys` as the `keytype` argument.

`keys` returns keys from the `Mart` of the type specified by it's
`keytype` argument.

`select` is meant to be used with these other methods and has arguments
that take the kinds of values that these other methods return. `select`
will retrieve the results as a data.frame based on parameters for
selected `keys` and `columns` and `keytype` arguments.

## Author

Marc Carlson

## Examples

``` r
if (FALSE) { # interactive()
## 1st create a Mart object and specify the dataset
mart <- useEnsembl(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)
## you can list the keytypes
keytypes(mart)
## you can list the columns
columns(mart)
## And you can extract keys when this is supported for your keytype of interest
k <- keys(mart, keytype="chromosome_name")
head(k)
## You can even do some pattern matching on the keys
k <- keys(mart, keytype="chromosome_name", pattern="LRG")
head(k)
## Finally you can use select to extract records for things that you are
## interested in.
affy <- c("202763_at", "209310_s_at", "207500_at")
select(mart, keys=affy, columns=c('affy_hg_u133_plus_2','entrezgene_id'),
  keytype='affy_hg_u133_plus_2')
}
```
