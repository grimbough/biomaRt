test_that("exportFASTA() on gene-sequence data.frame", {
  test_df <- data.frame(
    sequence = c(strrep("ATCG", 10), strrep("GCTA", 10)),
    gene = c("gene1", "gene2")
  )

  expect_snapshot(exportFASTA(test_df, stdout()))
})

test_that("exportFASTA() on 4 columns data", {
  test_df <- data.frame(
    chromosome = c("chr1", "chr1"),
    start = c(1, 55),
    end = c(1, 55) + 40,
    sequence = c(strrep("ATCG", 10), strrep("GCTA", 10))
  )

  expect_snapshot(exportFASTA(test_df, stdout()))
})
