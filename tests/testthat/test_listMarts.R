# library(biomaRt)
# context('Testing listMarts() function')
# 
# ## we expect this to error as the URL is out of date.
# expect_error(useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="www.biomart.org"), "Unable to determine the list of available marts")