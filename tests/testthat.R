library(testthat)
library(biomaRt)

#biomartCacheClear()
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = tempdir())

test_check("biomaRt", encoding = "UTF-8")