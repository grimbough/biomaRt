test_that("getGene() recognizes new ensembl name as ensembl", {
  # see https://github.com/Huber-group-EMBL/biomaRt/issues/122
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  expect_no_error(
    g <- getGene(id = "1939_at", type = "affy_hg_u95av2", mart = mart)
  )
  expect_identical(nrow(g), 1L)
})
