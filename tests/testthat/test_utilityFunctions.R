library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

ensembl = Mart(biomart = "ensembl", 
               dataset = "hsapiens_gene_ensembl",
               attributes = data.frame(
                   name = c("chromosome_name", "ensembl_gene_id"),
                   description = c("Chromosome/scaffold name", "Gene stable ID")
               ),
               filters = data.frame(
                   name = c("affy_hg_u133a_2", "chromosome_name", "transcript_tsl"),
                   description = c("AFFY HG U133A 2 probe ID(s) [e.g. 211600_at]", 
                                   "Chromosome/scaffold name",
                                   "Transcript Support Level (TSL)"),
                   type = c("id_list", "text", "boolean")
               )
               )

## create two test data.frames
bad_result <- data.frame("Not a real column" = 1:2, "Ensembl Gene ID" = 3:4, check.names = FALSE)
good_result <- data.frame("Chromosome/scaffold name" = 1:2, "Gene stable ID" = 3:4, check.names = FALSE)

test_that("Renaming columns - synthetic data", {
    ## this should warn that we can't match one of the column names (hopefully this never happens for real)
    expect_warning(.setResultColNames(result = bad_result, mart = ensembl, 
                                      attributes = c('chromosome_name', 'ensembl_gene_id')), 
                   "Problems assigning column names")
    ## check the reassignment of colnames from 'description' to 'name'
    expect_equal(colnames(.setResultColNames(result = good_result, mart = ensembl, 
                                             attributes = c('chromosome_name', 'ensembl_gene_id'))), 
                 c("chromosome_name", "ensembl_gene_id"))
    ## check we reorder them if needed
    expect_equal(colnames(.setResultColNames(result = good_result, mart = ensembl, 
                                             attributes = c('ensembl_gene_id', 'chromosome_name'))), 
                 c("ensembl_gene_id", "chromosome_name"))
})


#############################
context("Host name format processing")
#############################

test_that("URL formatting works", {
    ## adding http if needed
    host <- 'www.myurl.org'
    expect_equal(object = .cleanHostURL(host = host),
                 expected = "http://www.myurl.org")
    
    ## stripping trailing slash
    host <- 'http://www.myurl.org/'
    expect_equal(object = .cleanHostURL(host = host),
                 expected = "http://www.myurl.org")
    
    ## leave https already there
    host <- 'https://www.myurl.org'
    expect_equal(object = .cleanHostURL(host = host),
                 expected = "https://www.myurl.org")
    
    ## add 'www' to ensembl.org
    host <- 'ensembl.org'
    expect_equal(object = .cleanHostURL(host = host),
                 expected = "http://www.ensembl.org")
})

#############################
context("Query submission")
#############################

test_that("TSV and HTML result tables match", {
    host <- "https://www.ensembl.org:443/biomart/martservice?"
    query <- "<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query>
            <Query virtualSchemaName='default' uniqueRows='1' count='0' datasetConfigVersion='0.6' header='1' requestid='biomaRt' formatter='TSV'>
            <Dataset name = 'hsapiens_gene_ensembl'><Attribute name = 'ensembl_gene_id'/>
            <Attribute name = 'hgnc_symbol'/>
            <Filter name = \"ensembl_gene_id\" value = \"ENSG00000100036\" />
            </Dataset></Query>"
    expect_silent(res_tsv <- biomaRt:::.submitQueryXML(host, query))
    expect_silent(res_html <- biomaRt:::.fetchHTMLresults(host, query))
    expect_identical(res_html,
                     read.table(textConnection(res_tsv), sep = "\t", header = TRUE, 
                                check.names = FALSE, stringsAsFactors = FALSE))
})
