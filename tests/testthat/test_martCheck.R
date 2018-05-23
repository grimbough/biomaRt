library(biomaRt)
context('Testing martCheck function')

test_that("martCheck() catches bad input", { 
    expect_error(biomaRt:::martCheck())
    expect_error(biomaRt:::martCheck("INVALID_OBJECT"))
    
    ensembl <- useMart("ensembl")
    expect_error(biomaRt:::martCheck(ensembl), 
                 regex = "No dataset selected, please select a dataset first")
    
    ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
    expect_error(biomaRt:::martCheck(ensembl, biomart = "Not_real_mart"), 
                 regexp = "This function only works when used with")
})

test_that('martCheck() is quiet for valid input', {
    expect_silent(biomaRt:::martCheck(ensembl, 
                                      biomart = "ENSEMBL_MART_ENSEMBL"))
})
