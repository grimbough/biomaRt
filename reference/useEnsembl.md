# Connects to the selected BioMart database and dataset hosted by Ensembl

A first step in using the biomaRt package is to select a BioMart
database and dataset to use. The `useEnsembl()` function enables one to
connect to a specified BioMart database and dataset hosted by Ensembl
without having to specify the Ensembl URL. To know which BioMart
databases are available see the
[`listEnsembl()`](https://huber-group-embl.github.io/biomaRt/reference/listEnsembl.md)
and
[`listEnsemblGenomes()`](https://huber-group-embl.github.io/biomaRt/reference/listEnsembl.md)
functions. To know which datasets are available within a BioMart
database, first select the BioMart database using `useEnsembl()` and
then use the
[`listDatasets()`](https://huber-group-embl.github.io/biomaRt/reference/listDatasets.md)
function on the selected Mart object.

## Usage

``` r
useEnsembl(
  biomart,
  dataset,
  host,
  version = NULL,
  GRCh = NULL,
  mirror = NULL,
  verbose = FALSE
)

useEnsemblGenomes(biomart, dataset, host = NULL)
```

## Arguments

- biomart:

  BioMart database name you want to connect to. Possible database names
  can be retrieved with the function
  [`listEnsembl()`](https://huber-group-embl.github.io/biomaRt/reference/listEnsembl.md)

- dataset:

  Dataset you want to use. To see the different datasets available
  within a biomaRt you can e.g. do: mart = useEnsembl('genes'), followed
  by listDatasets(mart).

- host:

  Host to connect to. Only needs to be specified if different from
  www.ensembl.org. For `useEnsemblGenomes()` this argument can be used
  to specify an archive site.

- version:

  Ensembl version to connect to when wanting to connect to an archived
  Ensembl version

- GRCh:

  GRCh version to connect to if not the current GRCh38, currently this
  can only be 37

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  'www', 'useast', 'asia'. If no mirror is specified the primary site at
  www.ensembl.org will be used. Mirrors are not available for the
  Ensembl Genomes databases.

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

The `mirror` argument can be considered as a "preferred choice" when
connecting to Ensembl. If the argument is provided then connectivity to
that mirror will be tested. If it responds positively then the requested
mirror will be used. If the response is a failure each of the remaining
mirrors will be selected at random and tested until a working server is
found. Once identified that Ensembl server will be associated with the
returned `Mart` object and will be used for all queries.

## Author

Steffen Durinck & Mike Smith

## Examples

``` r
if (FALSE) { # interactive()
mart <- useEnsembl("ENSEMBL_MART_ENSEMBL")

## using the US East mirror
us_mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", mirror = "useast")

## using the Arabidopsis thaliana genes dataset in Ensembl Plants
plants_mart <- useEnsemblGenomes(
  biomart = "plants_mart",
  dataset = "athaliana_eg_gene"
)

## using the Cucumis melo genes dataset in the Ensembl Plants 56 archive
plants_mart <- useEnsemblGenomes(
  biomart = "plants_mart",
  dataset = "cmelo_eg_gene",
  host = "https://feb2023-plants.ensembl.org/"
)
}
```
