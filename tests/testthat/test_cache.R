library(biomaRt)

context("Result caching")

go <- c("GO:0051330","GO:0000080","GO:0000114","GO:0000082")
chrom <- c(17,20,"Y")
attributes <- "hgnc_symbol"
filters <- c("go","chromosome_name")
values <- list(go, chrom)
mart <- useEnsembl("ensembl", "hsapiens_gene_ensembl", mirror = "uswest")

test_that("Hashing is order insensitive", {
   expect_identical(
       biomaRt:::.createHash(mart, attributes, filters, values),
       biomaRt:::.createHash(mart, rev(attributes), rev(filters), rev(values))
   )
})
