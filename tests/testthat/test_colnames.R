library(biomaRt)

ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)

context('Testing column name assignments generation')
## create two test data.frames
bad_result <- data.frame("Not a real column" = 1:2, "Ensembl Gene ID" = 3:4, check.names = FALSE)
good_result <- data.frame("Chromosome/scaffold name" = 1:2, "Gene ID" = 3:4, check.names = FALSE)
## this should warn that we can't match one of the column names (hopefully this never happens for real)
expect_warning(.setResultColNames(result = bad_result, mart = ensembl, attributes = c('chromosome_name', 'ensembl_gene_id')), "Problems assigning column names")
## check the reassignment of colnames from 'description' to 'name'
expect_equal(colnames(.setResultColNames(result = good_result, mart = ensembl, attributes = c('chromosome_name', 'ensembl_gene_id'))), 
             c("chromosome_name", "ensembl_gene_id"))
## check we reorder them if needed
expect_equal(colnames(.setResultColNames(result = good_result, mart = ensembl, attributes = c('ensembl_gene_id', 'chromosome_name'))), 
             c("ensembl_gene_id", "chromosome_name"))

## check we can handle the ambiguous description field in some datasets
attributes=c("ensembl_transcript_id",
             "cdna", 
             ## Description is 'Query protein or transcript ID' lots of matches
             "meugenii_homolog_canonical_transcript_protein")
mart <- useMart(biomart = "ensembl",
              host = "www.ensembl.org",
              dataset ="mmusculus_gene_ensembl")
seq <- getBM(filter = "ensembl_gene_id",
           values = "ENSMUSG00000000103",
           attributes = attributes,
           mart = mart)
expect_equal(colnames(seq), attributes)