library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

plants <- Mart(biomart = "plants_mart",
              dataset = "athaliana_eg_gene",
              host = "https://plants.ensembl.org")

mouse <- Mart(biomart = "mouse_genes",
              dataset = "mmusculus_gene_ensembl",
              host = "https://www.ensembl.org")

mouse_strains <- Mart(biomart = "mouse_strains",
              dataset = "mmc57bl6nj_gene_ensembl",
              host = "https://www.ensembl.org")

test_that("Error with separate hosts", { 
    expect_error(getLDS(attributes="ensembl_gene_id", mart=mouse, 
                        attributesL="ensembl_gene_id", martL=plants),
                 regexp = 'Both datasets must be located on the same host')
})

test_that("We get an error with different Marts on the same host", { 
    expect_error(getLDS(attributes="ensembl_gene_id", mart=mouse, 
                        attributesL="ensembl_gene_id", martL=mouse_strains),
                 regexp = 'Both datasets must be located in the same Mart')
})

