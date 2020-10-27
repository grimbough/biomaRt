library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

context('Testing martCheck function')

test_that("martCheck() catches bad input", { 
    expect_error(biomaRt:::martCheck())
    expect_error(biomaRt:::martCheck("INVALID_OBJECT"))
    
    ensembl <- Mart(biomart = "ensembl")
    
    expect_error(biomaRt:::martCheck(ensembl), 
                 regex = "No dataset selected, please select a dataset first")
    
    ensembl <- Mart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    expect_error(biomaRt:::martCheck(ensembl, biomart = "Not_real_mart"), 
                 regexp = "This function only works when used with")
})

