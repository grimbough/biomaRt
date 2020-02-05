library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

context("useMart() functionality")

## checking the show() method
ensembl <- useMart("ensembl")
ensembl_with_dataset <- useDataset(ensembl, 
                                   dataset = "xtropicalis_gene_ensembl")

test_that("Show gives sensible dataset information", {
    expect_output(object = show(ensembl), 
                  regexp = "No dataset selected")
    expect_output(object = show(ensembl_with_dataset), 
                  regexp = "Using the xtropicalis_gene_ensembl dataset")
})



test_that("Deprecation warning produced", {
    
    expect_warning(useMart(biomart = "ensembl", host="www.ensembl.org", ensemblRedirect = FALSE),
                   regexp = "The argument \"ensemblRedirect\" has been deprecated and will be removed")
    
})
