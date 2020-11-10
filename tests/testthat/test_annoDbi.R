
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
                   type = c("id_list", "text", "boolean"),
                   options = c("[]", "[1,2,3,4,CHR_HG1_PATCH]", "[only,excluded]")
               )
)

test_that("AnnotationDbi style interface works", {
    
    expect_identical(keytypes(ensembl), martFilters(ensembl)$name)
    
    expect_identical(columns(ensembl), martAttributes(ensembl)$name)
    
    expect_identical(keys(ensembl, "chromosome_name"),
                     listFilterOptions(ensembl, filter = "chromosome_name"))
    
})