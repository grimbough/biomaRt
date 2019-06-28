library(biomaRt)

context("Result caching")

go <- c("GO:0051330","GO:0000080","GO:0000114","GO:0000082")
chrom <- c(17,20,"Y")
attributes <- "hgnc_symbol"
filters <- c("go","chromosome_name")
values <- list(go, chrom)
ensembl <- useEnsembl("ensembl", "hsapiens_gene_ensembl", mirror = "uswest")

test_that("Hashing is order insensitive", {
   expect_identical(
       biomaRt:::.createHash(ensembl, attributes, filters, values),
       biomaRt:::.createHash(ensembl, rev(attributes), rev(filters), rev(values))
   )
})

test_that("Cache details are printed", {
    expect_message( biomartCacheInfo(),
                   regexp = "biomaRt cache")  
})

test_that("Cache can be cleared", {
    cache_file <- rappdirs::user_cache_dir(appname = "biomaRt")
    expect_true( file.exists( cache_file ) )
    expect_silent( biomaRt:::biomartCacheClear() )
    expect_false( file.exists( cache_file) )
})