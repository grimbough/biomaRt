library(biomaRt)
context('getLDS() function')

ensembl_hsapiens <- useMart("ENSEMBL_MART_ENSEMBL", 
                dataset="hsapiens_gene_ensembl")
plants <- useMart("plants_mart", host="plants.ensembl.org", 
                 dataset="athaliana_eg_gene")
ensembl_rnorvegicus <- useMart("ENSEMBL_MART_ENSEMBL", 
                       dataset="rnorvegicus_gene_ensembl")

test_that("Error with separate hosts", { 
    expect_error(getLDS(attributes="ensembl_gene_id", mart=ensembl_hsapiens, 
                        attributesL="ensembl_gene_id", martL=plants),
                 regexp = 'Both datasets must be located on the same host')
})


test_that("Find human/rat homologs", { 
    
    expect_is(x <- getLDS(attributes="ensembl_gene_id", mart=ensembl_hsapiens, 
                          filters = "ensembl_gene_id", values = "ENSG00000084453",
                          attributesL="ensembl_gene_id", martL=ensembl_rnorvegicus),
              "data.frame")
    expect_true("ENSRNOG00000047493" %in% x[,2])
})
