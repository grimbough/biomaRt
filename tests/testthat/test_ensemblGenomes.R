test_that("Error handling works", {
  
    expect_error(useEnsemblGenomes(), 
                 regex = "You must provide the argument 'biomart'")
    
    with_mock(
        listEnsemblGenomes = function(includeHosts) { 
            data.frame(
                biomart = c("protists_mart", "fungi_mart"),
                version = c("Ensembl Protists Genes 48", "Ensembl Fungi Genes 48")
                )
            },
    expect_error(useEnsemblGenomes(biomart = "NOT_A_REAL_MART"),
                 regexp = "is not in the list of available Marts")
    )
    expect_true(is(useEnsemblGenomes("protists_mart", "bnatans_eg_gene"), "Mart"))
      
})