library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

test_that("useMart can connect to ensembl", {
  expect_silent(ensembl <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl'))
})

#############################
## checking the show() method
#############################

test_that("show() reports missing dataset", {
        expect_output(object = show( Mart(biomart = "ensembl") ), 
                      regexp = "No dataset selected")
})

###############

test_that("show() reports dataset name if present", {
    expect_output(object = show( Mart(biomart = "ensembl", dataset = "xtropicalis_gene_ensembl" )), 
                  regexp = "Using the xtropicalis_gene_ensembl dataset")
})

#############

test_that("Deprecation warning produced", {
    expect_error(useMart(biomart = "ensembl", host="www.ensembl.org", ensemblRedirect = FALSE))
})



