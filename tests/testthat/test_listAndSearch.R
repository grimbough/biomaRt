library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

context("list and search functions")

ensembl <- useMart("ensembl")
ensembl_with_dataset <- useDataset(ensembl, 
                                   dataset = "xtropicalis_gene_ensembl")

test_that("Fail with no mart argument", {
    expect_error(searchDatasets(), "Argument 'mart' must be specified")
    expect_error(searchAttributes(), "Argument 'mart' must be specified")
    expect_error(searchFilters(), "Argument 'mart' must be specified")
})

test_that("'Long' table of results for no search term", {
    
    expect_is(x <- searchDatasets(ensembl), class = 'data.frame')
    expect_gt(nrow(x), 20) ## arbitrary definition of 'long'
    
    expect_is(x <- searchAttributes(ensembl_with_dataset), class = 'data.frame')
    expect_gt(nrow(x), 20) 
    
    expect_is(x <- searchFilters(ensembl_with_dataset), class = 'data.frame')
    expect_gt(nrow(x), 20) 
})

test_that("Message when nothing found", {
    
    expect_message(searchDatasets(ensembl, pattern = "foobaa"), "No matching") %>%
        expect_null()

    expect_message(searchAttributes(ensembl_with_dataset, pattern = "foobaa"), "No matching") %>%
        expect_null()
    
    expect_message(searchFilters(ensembl_with_dataset, pattern = "foobaa"), "No matching") %>%
        expect_null()
})


test_that("'Sensible' table of results for specific search term", {
    
    expect_is(x <- searchDatasets(ensembl, pattern = "hsapiens"), class = 'data.frame')
    expect_equal(nrow(x), 1) ## only one human dataset
    
    expect_is(x <- searchAttributes(ensembl_with_dataset, pattern = "ensembl.*id$"), 
              class = 'data.frame')
    expect_gt(nrow(x), 0) 
    expect_lt(nrow(x), 25) ## arbitrary definition of 'too long'
    
    expect_is(x <- searchFilters(ensembl_with_dataset, pattern = "ensembl.*id$"), 
              class = 'data.frame')
    expect_gt(nrow(x), 0) 
    expect_lt(nrow(x), 10) ## arbitrary definition of 'too long'
})