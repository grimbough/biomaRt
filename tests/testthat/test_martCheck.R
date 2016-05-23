library(biomaRt)
context('Testing martCheck function')

expect_error(biomaRt:::martCheck())
expect_error(biomaRt:::martCheck("INVALID_OBJECT"))

ensembl=useMart("ensembl")
expect_error(biomaRt:::martCheck(ensembl), "No dataset selected, please select a dataset first")

ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
expect_error(biomaRt:::martCheck(ensembl, biomart = "Not_real_mart"), "This function only works when used with the Not_real_mart BioMart")
expect_silent(biomaRt:::martCheck(ensembl, biomart = "ENSEMBL_MART_ENSEMBL"))
