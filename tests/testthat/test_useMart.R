library(biomaRt)

## checking the show() method
ensembl <- useMart("ensembl")
ensembl_with_dataset <- useDataset(ensembl, 
                                   dataset = "xtropicalis_gene_ensembl")

test_that("Show give sensible dataset information", {
    expect_output(object = show(ensembl), 
                  regexp = "No dataset selected")
    expect_output(object = show(ensembl_with_dataset), 
                  regexp = "Using the xtropicalis_gene_ensembl dataset")
})
