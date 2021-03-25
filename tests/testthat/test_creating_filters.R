ensembl <- Mart(biomart = "ensembl", 
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

test_that("Filter generation works", {
    expect_equal(biomaRt:::.generateFilterXML(filters = c('affy_hg_u133a_2', 'chromosome_name'),
                                              values = list(affyid=c('1939_at','1000_at'), chromosome= '16'),
                                              mart = ensembl)[[1]],
                 "<Filter name = \"affy_hg_u133a_2\" value = \"1939_at,1000_at\" /><Filter name = \"chromosome_name\" value = \"16\" />")
    
    expect_equal(biomaRt:::.generateFilterXML(filters = 'chromosome_name',
                                              values = '16',
                                              mart = ensembl)[[1]],
                 "<Filter name = \"chromosome_name\" value = \"16\" />")
    
    expect_equal(biomaRt:::.generateFilterXML(filters = 'chromosome_name',
                                              values = c('16', '18'),
                                              mart = ensembl)[[1]],
                 "<Filter name = \"chromosome_name\" value = \"16,18\" />")
    
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
                 "<Filter name = \"transcript_tsl\" excluded = \"0\" />")
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
                 "<Filter name = \"chromosome_name\" value = \"16,18\" />")
})

test_that("numeric values to filters work", {
    expect_equal(biomaRt:::.generateFilterXML(filters = 'chromosome_name',
                                              values = c(16, 18),
                                              mart = ensembl)[[1]],
                 "<Filter name = \"chromosome_name\" value = \"16,18\" />")
})

# test_that("large numbers of values are split into chunks", {
#     
#     for(max_size in c(2, 10, 50, 100, 200, 500)) {
#         expect_length(biomaRt:::.generateFilterXML(filters = 'start',
#                                  values = 1:500, maxChunkSize = max_size,
#                                  mart = ensembl),
#                       ceiling(500 / max_size))
#     }
# })
