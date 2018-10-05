library(biomaRt)

ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)

#############################
context('Column name assignments')
#############################

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

test_that("Renaming columns - real data", {        
    ## check we can handle the ambiguous description field in some datasets
    attributes=c("ensembl_transcript_id",
                 ## Description is 'Query protein or transcript ID' lots of matches
                 "neugenii_homolog_canonical_transcript_protein")
    
    mart <- useMart(biomart = "ensembl",
                    host = "https://www.ensembl.org",
                    dataset ="mmusculus_gene_ensembl",
                    port = 443)
    
    res <- getBM(filter = "ensembl_gene_id",
                 values = "ENSMUSG00000028798",
                 attributes = attributes,
                 mart = mart)
    expect_equal(colnames(res), attributes)
    
})

#############################
context('Testing filter XML generation')
#############################

test_that("Filter generation works", {
    expect_equal(biomaRt:::.generateFilterXML(filters = c('affy_hg_u133a_2', 'chromosome_name'),
                                    values = list(affyid=c('1939_at','1000_at'), chromosome= '16'),
                                    mart = ensembl)[[1]],
                 "<Filter name = 'affy_hg_u133a_2' value = '1939_at,1000_at' /><Filter name = 'chromosome_name' value = '16' />")
    
    expect_equal(biomaRt:::.generateFilterXML(filters = 'chromosome_name',
                                    values = '16',
                                    mart = ensembl)[[1]],
                 "<Filter name = 'chromosome_name' value = '16' />")
    
    expect_equal(biomaRt:::.generateFilterXML(filters = 'chromosome_name',
                                    values = c('16', '18'),
                                    mart = ensembl)[[1]],
                 "<Filter name = 'chromosome_name' value = '16,18' />")
    
    expect_equal(biomaRt:::.generateFilterXML(filters = ''), "")
})

test_that("Boolean filter handled correctly", {
    expect_error(biomaRt:::.generateFilterXML(filters = 'transcript_tsl',
                                values = '16',
                                mart = ensembl),
             "'transcript_tsl' is a boolean filter")
    expect_equal(biomaRt:::.generateFilterXML(filters = 'transcript_tsl',
                                values = TRUE,
                                mart = ensembl)[[1]],
             "<Filter name = 'transcript_tsl' excluded = \"0\" />")
})

test_that("list is required for multiple filters", {
    
    expect_error(.generateFilterXML(filters = c('affy_hg_u133a_2', 'chromosome_name'),
                       values = c(affyid=c('1939_at','1000_at'), chromosome= '16'),
                       mart = ensembl)[[1]])
})

test_that("passing a single column data.frame works", {
    expect_equal(biomaRt:::.generateFilterXML(filters = 'chromosome_name',
                                              values = data.frame('chr' = c('16', '18')),
                                              mart = ensembl)[[1]],
                 "<Filter name = 'chromosome_name' value = '16,18' />")
})

test_that("numeric values to filters work", {
    expect_equal(biomaRt:::.generateFilterXML(filters = 'chromosome_name',
                                              values = c(16, 18),
                                              mart = ensembl)[[1]],
                 "<Filter name = 'chromosome_name' value = '16,18' />")
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