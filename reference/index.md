# Package index

## `Mart` class and methods

- [`show(`*`<Mart>`*`)`](https://huber-group-embl.github.io/biomaRt/reference/Mart-class.md)
  : Class Mart
- [`keys(`*`<Mart>`*`)`](https://huber-group-embl.github.io/biomaRt/reference/select-methods.md)
  [`keytypes(`*`<Mart>`*`)`](https://huber-group-embl.github.io/biomaRt/reference/select-methods.md)
  [`columns(`*`<Mart>`*`)`](https://huber-group-embl.github.io/biomaRt/reference/select-methods.md)
  [`select(`*`<Mart>`*`)`](https://huber-group-embl.github.io/biomaRt/reference/select-methods.md)
  : Retrieve information from the BioMart databases

## Create or update `Mart` objects

- [`listMarts()`](https://huber-group-embl.github.io/biomaRt/reference/listMarts.md)
  : lists the available BioMart databases
- [`listEnsembl()`](https://huber-group-embl.github.io/biomaRt/reference/listEnsembl.md)
  [`listEnsemblGenomes()`](https://huber-group-embl.github.io/biomaRt/reference/listEnsembl.md)
  : lists the available BioMart databases hosted by Ensembl
- [`listEnsemblArchives()`](https://huber-group-embl.github.io/biomaRt/reference/listEnsemblArchives.md)
  : Lists the available archived versions of Ensembl
- [`useDataset()`](https://huber-group-embl.github.io/biomaRt/reference/useDataset.md)
  : Select a dataset to use and updates Mart object
- [`useEnsembl()`](https://huber-group-embl.github.io/biomaRt/reference/useEnsembl.md)
  [`useEnsemblGenomes()`](https://huber-group-embl.github.io/biomaRt/reference/useEnsembl.md)
  : Connects to the selected BioMart database and dataset hosted by
  Ensembl
- [`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md)
  : Connects to the selected BioMart database and dataset

## Extract data from BioMart databases

- [`listAttributes()`](https://huber-group-embl.github.io/biomaRt/reference/listAttributes.md)
  [`searchAttributes()`](https://huber-group-embl.github.io/biomaRt/reference/listAttributes.md)
  : lists the attributes available in the selected dataset
- [`attributePages()`](https://huber-group-embl.github.io/biomaRt/reference/attributePages.md)
  : Gives a summary of the attribute pages
- [`listDatasets()`](https://huber-group-embl.github.io/biomaRt/reference/listDatasets.md)
  [`searchDatasets()`](https://huber-group-embl.github.io/biomaRt/reference/listDatasets.md)
  : List or search the datasets available in the selected BioMart
  database
- [`listFilters()`](https://huber-group-embl.github.io/biomaRt/reference/listFilters.md)
  [`searchFilters()`](https://huber-group-embl.github.io/biomaRt/reference/listFilters.md)
  : List or search the filters available in the selected dataset
- [`searchFilterOptions()`](https://huber-group-embl.github.io/biomaRt/reference/listFilterOptions.md)
  [`listFilterOptions()`](https://huber-group-embl.github.io/biomaRt/reference/listFilterOptions.md)
  : List or search the options available for a specified filter.
- [`filterType()`](https://huber-group-embl.github.io/biomaRt/reference/filterType.md)
  : Displays the filter type
- [`getBM()`](https://huber-group-embl.github.io/biomaRt/reference/getBM.md)
  : Retrieves information from the BioMart database
- [`getGene()`](https://huber-group-embl.github.io/biomaRt/reference/getGene.md)
  : Retrieves gene annotation information given a vector of identifiers
- [`getLDS()`](https://huber-group-embl.github.io/biomaRt/reference/getLDS.md)
  : Retrieves information from two linked datasets
- [`getSequence()`](https://huber-group-embl.github.io/biomaRt/reference/getSequence.md)
  : Retrieves sequences

## Other utilities

These functions provide functionality previously covered by BioMart, but
no longer available in recent Ensembl releases.

- [`getHomologs()`](https://huber-group-embl.github.io/biomaRt/reference/getHomologs.md)
  : List homologous genes between two species.

## Data export

- [`exportFASTA()`](https://huber-group-embl.github.io/biomaRt/reference/exportFASTA.md)
  : Exports getSequence results to FASTA format

## Network and cache utilities

- [`setEnsemblSSL()`](https://huber-group-embl.github.io/biomaRt/reference/setEnsemblSSL.md)
  : Save system specific SSL settings for contacting Ensembl

- [`biomartCacheClear()`](https://huber-group-embl.github.io/biomaRt/reference/biomartCache.md)
  [`biomartCacheInfo()`](https://huber-group-embl.github.io/biomaRt/reference/biomartCache.md)
  :

  biomaRt result caching

## Dataset

- [`ensembl_versions`](https://huber-group-embl.github.io/biomaRt/reference/ensembl_versions.md)
  : Mapping of Ensembl versions to their release dates
- [`NP2009code()`](https://huber-group-embl.github.io/biomaRt/reference/NP2009code.md)
  : Display the analysis code from the 2009 Nature protocols paper

## Deprecated of Defunct functions

- [`biomaRt-deprecated`](https://huber-group-embl.github.io/biomaRt/reference/biomaRt-deprecated.md)
  [`filterOptions`](https://huber-group-embl.github.io/biomaRt/reference/biomaRt-deprecated.md)
  [`searchFilterValues`](https://huber-group-embl.github.io/biomaRt/reference/biomaRt-deprecated.md)
  [`listFilterValues`](https://huber-group-embl.github.io/biomaRt/reference/biomaRt-deprecated.md)
  :

  Deprecated and defunct functions in package biomaRt
