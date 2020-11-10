cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

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

test_that("getBM() gives appropriate error messages", {
    
    expect_error(getBM(), 
                 "You must provide a valid Mart object")
    
    mart_with_no_dataset <- Mart(biomart = "ensembl")
    expect_error(getBM(mart = mart_with_no_dataset), 
                 "No dataset selected, please select a dataset first")
    
    expect_error(getBM(mart = ensembl), 
                 "Argument 'attributes' must be specified")
    
    expect_error(getBM(mart = ensembl,
                       attributes = "ensembl_gene_id",
                       filters = list("chromosome_name")), 
                 "Argument 'filters' must be a named list when sent as a list")
    

})


test_that("getBMlist removal message is shown", {
    expect_error(getBMlist(),
                 class = "defunctError")
})



test_that("getBM returns sensible things", {
    
    stub(getBM,
         '.submitQueryXML',
         'Gene stable ID\tChromosome/scaffold name\nENSG01\t13\nENSG02\t15\nENSG03\t17\n'
    )
    expect_is(getBM(mart = ensembl,
                    attributes = c("ensembl_gene_id", "chromosome_name")),
              "data.frame")
    
    expect_warning(getBM(mart = ensembl,
                         attributes = c("ensembl_gene_id", "chromosome_name"),
                         filters = list(chromosome_name = c(1)),
                         values = "1"), 
                   "Argument 'values' should not be used when argument 'filters' is a list and will be ignored")
    
})

#######################
## the definition_1006 entry for this gene includes unescaped new lines, so the HTML result is requested
########################
# test_that("HTML reading code is used when needed", {
#     expect_silent(ensembl <- useEnsembl("ensembl", 
#                                         dataset = 'hsapiens_gene_ensembl',
#                                         mirror = "www"))
#     attributes <-  c("ensembl_gene_id", "go_id", "definition_1006")
#     expect_silent(go_sets <- getBM(attributes =  attributes,
#                                    filters = "ensembl_gene_id",
#                                    values = c('ENSG00000100036'),
#                                    mart = ensembl,
#                                    bmHeader = FALSE,
#                                    useCache = FALSE))
#     expect_is(go_sets, "data.frame")
#     expect_equal(colnames(go_sets), attributes)
# })
 
