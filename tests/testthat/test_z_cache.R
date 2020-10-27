library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

context("Result caching")

go <- c("GO:0051330","GO:0000080","GO:0000114","GO:0000082")
chrom <- c(17,20,"Y")
attributes <- "hgnc_symbol"
filters <- c("go","chromosome_name")
values <- list(go, chrom)
ensembl <- Mart(biomart = "ensembl", 
                dataset = "hsapiens_gene_ensembl", 
                host = "www.ensembl.org")

test_that("Hashing is order insensitive", {
   expect_identical(
       biomaRt:::.createHash(ensembl, attributes, filters, values),
       biomaRt:::.createHash(ensembl, rev(attributes), rev(filters), rev(values))
   )
})

test_that("Environment variable for cache location is used", {
    expect_message(biomartCacheInfo(),
                   regexp = "biomart_cache_test")
})



## We should refactor the code for adding entries to the cache, then modify
## this test to check it works correctly.
# test_that("We find cache for previous query", {
#     
#     ## construct the same settings as a query used in test_utilityFunctions.R
#     mart <- Mart(
#         biomart = "ENSEMBL_MART_ENSEMBL",
#         dataset = "mmusculus_gene_ensembl",
#         host = "https://useast.ensembl.org:443/biomart/martservice?redirect=no")
#     filters <- "ensembl_gene_id"
#     values <- "ENSMUSG00000028798"
#     attributes <- c("ensembl_transcript_id", 
#                     "neugenii_homolog_canonical_transcript_protein")
#     hash <- biomaRt:::.createHash(mart, attributes, filters, values)
#     bfc <- BiocFileCache::BiocFileCache(biomartCacheInfo(), ask = FALSE)
#     
#     expect_true(biomaRt:::.checkCache(bfc, hash))
# 
# })

test_that("Cache details are printed", {
    expect_message( biomartCacheInfo(),
                   regexp = "biomaRt cache")  
})

# test_that("Cache can be cleared", {
#     cache_file <- biomartCacheInfo()
#     expect_true( file.exists( cache_file ) )
#     expect_silent( biomartCacheClear() )
#     expect_false( file.exists( cache_file) )
# })