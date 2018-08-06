library(biomaRt)

########################
context('getBM()')
########################

test_that("Fail with no arguments", {
    expect_error(getBM(), "You must provide a valid Mart object")
})

test_that("Fail when no dataset is specified", {
     ensembl=useMart("ensembl")
     expect_error(getBM(mart = ensembl), "No dataset selected, please select a dataset first")
})
# 
# ## no attributes are specified')
# ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
# expect_error(getBM(mart = ensembl), "Argument 'attributes' must be specified")
# 
# t1 <- getBM(attributes='entrezgene', filters = 'affy_hg_u133_plus_2', values = '207500_at', mart = ensembl)[1,1]
# expect_equal(as.integer(t1), 838)
# 
# 
 
