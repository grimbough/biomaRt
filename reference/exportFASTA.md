# Exports getSequence results to FASTA format

Exports getSequence results to FASTA format

## Usage

``` r
exportFASTA(sequences, file)
```

## Arguments

- sequences:

  A data.frame that was the output of the
  [`getSequence()`](https://huber-group-embl.github.io/biomaRt/reference/getSequence.md)
  function

- file:

  File to which you want to write the data

## Author

Steffen Durinck

Hugo Gruson

## Examples

``` r
if (FALSE) { # interactive()
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

seq <- getSequence(
  id = "BRCA1",
  type = "hgnc_symbol",
  seqType = "cdna",
  mart = mart
)
exportFASTA(seq, file = "test.fasta")
}
```
