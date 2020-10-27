library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

########################
context('getBM()')
########################

test_that("Fail with no arguments", {
    expect_error(getBM(), 
                 "You must provide a valid Mart object")
})

test_that("Fail when no dataset is specified", {
     mart_with_no_dataset <- Mart(biomart = "ensembl")
     expect_error(getBM(mart = mart_with_no_dataset), 
                  "No dataset selected, please select a dataset first")
})


#######################
## the definition_1006 entry for this gene includes unescaped new lines, so the HTML result is requested
########################
test_that("HTML reading code is used when needed", {
    expect_silent(ensembl <- useEnsembl("ensembl", 
                                        dataset = 'hsapiens_gene_ensembl',
                                        mirror = "www"))
    attributes <-  c("ensembl_gene_id", "go_id", "definition_1006")
    expect_silent(go_sets <- getBM(attributes =  attributes,
                                   filters = "ensembl_gene_id",
                                   values = c('ENSG00000100036'),
                                   mart = ensembl,
                                   bmHeader = FALSE,
                                   useCache = FALSE))
    expect_is(go_sets, "data.frame")
    expect_equal(colnames(go_sets), attributes)
})
 
