library(biomaRt)
context('Testing getBM() function')

## Test: missing mart')
expect_error(getBM(), "You must provide a valid Mart object")

## Test: no dataset is specified')
ensembl=useMart("ensembl")
expect_error(getBM(mart = ensembl), "No dataset selected, please select a dataset first")

## no attributes are specified')
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
expect_error(getBM(mart = ensembl), "Argument 'attributes' must be specified")

expect_equal(getBM(attributes='entrezgene', filters = 'affy_hg_u133_plus_2', values = '207500_at', mart = ensembl)[1,1], 838)




context('Testing filter XML generation')
expect_equal(.generateFilterXML(filters = c('affy_hg_u133a_2', 'chromosome_name'), 
                   values = list(affyid=c('1939_at','1000_at'), chromosome= '16'), 
                   mart = ensembl),
             "<Filter name = 'affy_hg_u133a_2' value = '1939_at,1000_at' /><Filter name = 'chromosome_name' value = '16' />")

expect_equal(.generateFilterXML(filters = 'chromosome_name', 
                                values = '16', 
                                mart = ensembl),
             "<Filter name = 'chromosome_name' value = '16' />")

expect_equal(.generateFilterXML(filters = 'chromosome_name', 
                               values = c('16', '18'), 
                                mart = ensembl), 
             "<Filter name = 'chromosome_name' value = '16,18' />")

expect_equal(.generateFilterXML(filters = ''), "")

## testing with a boolean filter
expect_error(.generateFilterXML(filters = 'transcript_tsl',
                                values = '16',
                                mart = ensembl),
             "biomaRt error: transcript_tsl is a boolean filter")
expect_equal(.generateFilterXML(filters = 'transcript_tsl',
                                values = TRUE,
                                mart = ensembl),
             "<Filter name = 'transcript_tsl' excluded = \"0\" />")


context('Testing column name assignments generation')
bad_result <- data.frame("Not a real column" = 1:2, "Ensembl Gene ID" = 3:4, check.names = FALSE)
good_result <- data.frame("Chromosome Name" = 1:2, "Ensembl Gene ID" = 3:4, check.names = FALSE)
expect_warning(.setResultColNames(result = bad_result, mart = ensembl), "Problems assigning column names")
expect_equal(colnames(.setResultColNames(result = good_result, mart = ensembl)), 
             c("chromosome_name", "ensembl_gene_id"))


