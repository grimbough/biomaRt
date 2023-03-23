library(biomaRt)

cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

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
       .createHash(ensembl, attributes, filters, values),
       .createHash(ensembl, rev(attributes), rev(filters), rev(values))
   )
})

test_that("Environment variable for cache location is used", {
    expect_message(biomartCacheInfo(),
                   regexp = "biomart_cache_test")
})

## create an example hash and dataframe to test code with
hash <- .createHash(ensembl, attributes, filters, values)
result <- data.frame(
    name = c("affy_hg_u133a_2", "chromosome_name", "transcript_tsl"),
    description = c("AFFY HG U133A 2 probe ID(s) [e.g. 211600_at]", 
                    "Chromosome/scaffold name",
                    "Transcript Support Level (TSL)"),
    type = c("id_list", "text", "boolean"),
    options = c("[]", "[1,2,3,4,CHR_HG1_PATCH]", "[only,excluded]")
)
bfc <- BiocFileCache::BiocFileCache(biomartCacheInfo(), ask = FALSE)

test_that("Entries can be added to the cache", {
    expect_true(.addToCache(bfc = bfc, result = result, hash = hash))
})

test_that("We find cache for previous query", {
    expect_true(.checkInCache(bfc, hash = hash))
    expect_identical(result, .readFromCache(bfc, hash = hash))
})

## add an invalid file.  We are only expecting .rds files, and this is a 
## text file
tf <- tempfile()
writeLines(LETTERS, con = tf)
BiocFileCache::bfcadd(bfc, rname = "invalid-file", fpath = tf)
prev_count <- BiocFileCache::bfccount(bfc)

test_that("invalid cache entry detected", {
  expect_false(.checkValidCache(bfc, hash = "invalid-file"))
  expect_equal(BiocFileCache::bfccount(bfc), prev_count - 1)
})


## check we can't write multiple entries with the same hash
test_that("we can't create multiple entries with the same hash", {
  expect_true(.addToCache(bfc = bfc, result = 1:100, hash = "testing"))
  expect_false(.addToCache(bfc = bfc, result = 1:100, hash = "testing"))
})



test_that("Cache details are printed", {
    expect_message( biomartCacheInfo(),
                   regexp = "biomaRt cache")  
})

test_that("Cache can be cleared", {
    cache_file <- biomartCacheInfo()
    expect_true( file.exists( cache_file ) )
    expect_silent( biomartCacheClear() )
    expect_false( file.exists( cache_file) )
})