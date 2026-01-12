# List homologous genes between two species.

This function simplifies the querying of the Ensembl BioMart if you're
trying to return the homologs for one or more gene IDs between two
species.

## Usage

``` r
getHomologs(ensembl_gene_ids, species_from, species_to)
```

## Arguments

- ensembl_gene_ids:

  Character vector. This contains the Ensembl Gene IDs that you want to
  find the homologs for.

- species_from, species_to:

  Character vectors of length 1. These arguments specify the species the
  input IDs belong to (`species_from`) and the species you want to find
  the homologs in (`species_to`). These can be Ensembl genomes names
  e.g. "homo_sapiens" or "canis_lupus_familiaris" or common names e.g.
  "human" or "dog". The function will do it's best to parse common
  names, and will report and error if no match to an Ensembl genome can
  be made.

## Author

Mike Smith
