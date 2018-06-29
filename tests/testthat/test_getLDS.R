library(biomaRt)
context('getLDS() function')

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                dataset="hsapiens_gene_ensembl")
plants <- useMart("plants_mart", host="plants.ensembl.org", 
                 dataset="athaliana_eg_gene")

test_that("Error with separate hosts", { 
    expect_error(getLDS(attributes="ensembl_gene_id", mart=ensembl, 
                        attributesL="ensembl_gene_id", martL=plants),
                 regexp = 'Both datasets must be located on the same host')
})

