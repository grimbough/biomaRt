library(biomaRt)

###############################
context('listMarts()')
###############################

test_that("listMarts retruns a data.frame", {
    
    ensembl_marts <- listMarts(host = "www.ensembl.org")
    expect_is(ensembl_marts, class = "data.frame")
    expect_identical(colnames(ensembl_marts),
                     c("biomart", "version"))
    
})

test_that("Error when archive = TRUE", {
    
    expect_error(listMarts(host = "www.ensembl.org", archive = TRUE),
                 regexp = "Use listEnsemblArchives")
       
})

test_that("Error when old URL is used", {
    
    expect_error(listMarts(host="www.biomart.org"), 
                 "Unexpected format to the list of available marts")
    
})