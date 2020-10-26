library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

context("list and search functions")

ensembl <- useMart("ensembl")

## Construct example mart object without contacting Ensembl
ensembl_with_dataset <- Mart(biomart = "ensembl", 
                             dataset = "xtropicalis_gene_ensembl",
                             attributes = data.frame(name = c("A", "B", "C"),
                                                     description = c("AttrA", "AttrB", "AttrC"),
                                                     page = "P1"),
                             filters = data.frame(name = c("X", "Y", "Z"),
                                                  description = c("Filt_X", "Filt_Y", "Filt_Z")
                                 
                             ))

test_that("Fail with no mart argument", {
    expect_error(searchDatasets(), "Argument 'mart' must be specified")
    expect_error(searchAttributes(), "Argument 'mart' must be specified")
    expect_error(searchFilters(), "Argument 'mart' must be specified")
})

test_that("'Long' table of results for no search term", {
    
    expect_is(x <- searchDatasets(ensembl), class = 'data.frame')
    expect_gt(nrow(x), 20) ## arbitrary definition of 'long'
    
})

test_that("Return complete table of results for no search term", {
    
    expect_is(x <- searchAttributes(ensembl_with_dataset), class = 'data.frame')
    expect_identical(x, martAttributes(ensembl_with_dataset))  
    
    expect_is(x <- searchFilters(ensembl_with_dataset), class = 'data.frame')
    expect_identical(x, martFilters(ensembl_with_dataset)) 
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
    
    expect_is(x <- searchAttributes(ensembl_with_dataset, pattern = "Attr(A|B)"), 
              class = 'data.frame')
    expect_equal(nrow(x), 2) 

    expect_is(x <- searchFilters(ensembl_with_dataset, pattern = "Filt_Z"), 
              class = 'data.frame')
    expect_equal(nrow(x), 1) 
})