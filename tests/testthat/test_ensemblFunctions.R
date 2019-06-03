library(biomaRt)
context('Testing Ensembl specific functions')

test_that("useEnsembl() works", { 

  ensembl <- useEnsembl(biomart = "snp")
  expect_is(ensembl, class = "Mart")
  

  
})

test_that("useEnsembl() options", { 
  
  expect_silent(ensembl_asia <- useEnsembl(biomart = "snp", mirror = "uswest"))
  expect_equal(ensembl_asia@host, 
               "https://uswest.ensembl.org:443/biomart/martservice")
  
  expect_silent(ensembl_archive <- useEnsembl(biomart = "ensembl", version = 93))
  expect_equal(ensembl_archive@host,
               "http://jul2018.archive.ensembl.org:80/biomart/martservice")
})


test_that("useEnsembl() error handling", { 
  
  expect_error(useEnsembl(), regexp = "You must provide the argument")
  
  expect_warning(ensembl <- useEnsembl(biomart = "snp", mirror = "INVALID_MIRROR"), 
                 regexp = "Invalid mirror\\. Select a mirror")
  expect_equal(ensembl@host, 
               "https://www.ensembl.org:443/biomart/martservice")
  
  expect_error(useEnsembl(biomart = "snp", version = "00"), 
               regexp = "Specified Ensembl version is not available")
})
