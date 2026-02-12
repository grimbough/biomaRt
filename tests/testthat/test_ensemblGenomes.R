test_that("Error handling works", {
  skip_if_not_installed("mockery")

  expect_error(
    useEnsemblGenomes(),
    regexp = "You must provide the argument 'biomart'"
  )

  mockery::stub(
    useEnsemblGenomes,
    "listEnsemblGenomes",
    function(includeHosts, host) {
      data.frame(
        biomart = c("protists_mart", "fungi_mart"),
        version = c("Ensembl Protists Genes 48", "Ensembl Fungi Genes 48")
      )
    }
  )
  expect_error(
    useEnsemblGenomes(biomart = "plants_mart"),
    regexp = "is not in the list of available Marts"
  )
})

test_that("return type is correct", {
  skip_if_offline("protists.ensembl.org")
  skip_if_offline("fungi.ensembl.org")
  skip_if_offline("metazoa.ensembl.org")
  skip_if_offline("plants.ensembl.org")

  expect_s4_class(useEnsemblGenomes("protists_mart", "bnatans_eg_gene"), "Mart")
})
