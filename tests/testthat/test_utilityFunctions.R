cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)
bfc <- BiocFileCache::BiocFileCache(biomartCacheInfo(), ask = FALSE)

httr_config <- do.call(c, biomaRt:::.getEnsemblSSL())

ensembl <- Mart(biomart = "ensembl", 
               dataset = "hsapiens_gene_ensembl",
               host = "www.ensembl.org",
               attributes = data.frame(
                   name = c("chromosome_name", "ensembl_gene_id"),
                   description = c("Chromosome/scaffold name", "Gene stable ID")
               ),
               filters = data.frame(
                   name = c("affy_hg_u133a_2", "chromosome_name", "transcript_tsl"),
                   description = c("AFFY HG U133A 2 probe ID(s) [e.g. 211600_at]", 
                                   "Chromosome/scaffold name",
                                   "Transcript Support Level (TSL)"),
                   type = c("id_list", "text", "boolean"),
                   options = c("[]", "[1,2,3,4,CHR_HG1_PATCH]", "[only,excluded]")
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

example_return <- "GO term accession\tGene name\tGO domain\nGO:0004465\tLpl\tmolecular_function\nGO:0006631\tLpl\tbiological_process\nGO:0005509\tLpl\tmolecular_function\n"

test_that("Results processing works", {
    
    expect_is(result_table <- .processResults(postRes = example_return, mart = ensembl, 
                              fullXmlQuery = "", numAttributes = 3), 
              "data.frame")
    expect_identical(ncol(result_table), 3L)
    
    expect_error(.processResults(postRes = LETTERS, mart = ensembl, 
                                 fullXmlQuery = "", numAttributes = 3))
    expect_error(.processResults(postRes = example_return, mart = ensembl, 
                                 fullXmlQuery = "", numAttributes = 1))
    expect_error(.processResults(postRes = "Query ERROR", mart = ensembl, 
                                 fullXmlQuery = "", numAttributes = 3))
})


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
    ## We expect a warning about using http
    host <- 'ensembl.org'
    expect_warning(object = .cleanHostURL(host = host)) |>
      expect_equal(expected = "http://www.ensembl.org")
})

test_that("TSV and HTML result tables match", {
  
  skip_if_not_installed("mockery")
  
  host <- "https://www.ensembl.org:443/biomart/martservice?redirect=no"
  query <- "<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query>
            <Query virtualSchemaName='default' uniqueRows='1' count='0' datasetConfigVersion='0.6' header='1' requestid='biomaRt' formatter='TSV'>
            <Dataset name = 'hsapiens_gene_ensembl'><Attribute name = 'ensembl_gene_id'/>
            <Attribute name = 'hgnc_symbol'/>
            <Filter name = \"ensembl_gene_id\" value = \"ENSG00000100036\" />
            </Dataset></Query>"
  
  mockery::stub(where = .fetchHTMLresults, 
                what = ".submitQueryXML",
                how = '<?xml version="1.0"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
<head>
  <title></title>
  <link rel="stylesheet" type="text/css" href="/martview/martview.css" />
</head>
<body>

<table>
<tr>
  <th>Gene stable ID</th>
  <th>HGNC symbol</th>
</tr>
<tr>
  <td><a href="http://www.ensembl.org/homo_sapiens/Gene/Summary?db=core;g=ENSG00000100036" target="_blank">ENSG00000100036</a></td>
  <td><a href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:17058" target="_blank">SLC35E4</a></td>
</tr>

</table>
</body>
</html>')
  
  expect_silent(res_html <- .fetchHTMLresults(host, query, httr_config = httr_config))
  
  expect_identical(res_html, 
                  data.frame("Gene stable ID" = "ENSG00000100036", 
                             "HGNC symbol" = "SLC35E4",
                             check.names = FALSE)
                  )
})


test_that("http error codes are presented nicely", {
    expect_true(grepl(.createErrorMessage(error_code = 500, host = "test.com"),
                pattern = 'biomaRt has encountered an unexpected server error'))
    expect_true(grepl(.createErrorMessage(error_code = 509, host = "test.com"),
                pattern = 'biomaRt has exceeded the bandwidth allowance with this server'))
    
    expect_true(grepl(.createErrorMessage(error_code = 500, host = "www.ensembl.org")[2],
                      pattern = 'Consider trying one of the Ensembl mirrors'))
    expect_true(grepl(.createErrorMessage(error_code = 509, host = "www.ensembl.org")[2],
                      pattern = 'Consider trying one of the Ensembl mirrors'))
    
    ## testing an un-seen error code
    expect_true(grepl(.createErrorMessage(error_code = 999, host = "www.ensembl.org")[2],
                      pattern = 'Consider trying one of the Ensembl mirrors'))
    
})


test_that("we can search predefined filter values", {
    
    expect_equal(length(searchFilterOptions(ensembl, filter = "chromosome_name")), 5)
    
    expect_equal(length(searchFilterOptions(ensembl, "chromosome_name", "^[0-9]*$")), 4)
    
    expect_equal(searchFilterOptions(ensembl, "chromosome_name", "PATCH"), "CHR_HG1_PATCH")
})


test_that("deprecated functions show warnings", {
    
    expect_warning(searchFilterValues(ensembl, filter = "chromosome_name"))
    expect_warning(listFilterValues(ensembl, filter = "chromosome_name"))
})

test_that("attribute and filter tables are parsed correctly", {

    skip_if_not_installed('mockery')
  
    mockery::stub(.getAttrFilt, 
        'bmRequest',
        'ensembl_gene_id\tGene stable ID\tStable ID of the Gene\tfeature_page\thtml,txt,csv,tsv,xls\thsapiens_gene_ensembl__gene__main\tstable_id_1023\n',
    )
    expect_is(.getAttrFilt(mart = ensembl, verbose = TRUE, type = "attributes"), "data.frame")
    
    mockery::stub(.getAttributes, 
         '.getAttrFilt',
         read.table(text = "ensembl_gene_id\tGene stable ID\tStable ID of the Gene\tfeature_page\n",
                    sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
    )
    expect_is(.getAttributes(mart = ensembl, verbose = TRUE), "data.frame")
    
    
    mockery::stub(.getFilters, 
         '.getAttrFilt',
         read.table(text = "chromosome_name\tChromosome/scaffold name\t[]\t\tfilters\ttext\t=\tbnatans_eg_gene__gene__main\tname_1059\n\n",
                    sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
    )
    expect_is(.getFilters(mart = ensembl, verbose = TRUE), "data.frame")
      
})
