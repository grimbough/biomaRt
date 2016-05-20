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
expect_equal(.generateFilterXML(filters = c('affy_hg_u133a_2', 'chromosome'), 
                   values = list(affyid=c('1939_at','1000_at'), chromosome= '16'), 
                   mart = ensembl),
             "<Filter name = 'affy_hg_u133a_2' value = '1939_at,1000_at' /><Filter name = 'chromosome' value = '16' />")

expect_equal(.generateFilterXML(filters = 'chromosome', 
                                values = chromosome= '16', 
                                mart = ensembl),
             "<Filter name = 'chromosome' value = '16' />")