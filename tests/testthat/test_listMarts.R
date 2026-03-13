library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

with_mock_dir(
  "all_200",
  {
    test_that("listMarts returns a data.frame", {
      ensembl_marts <- listMarts(host = "https://www.ensembl.org")
      expect_s3_class(ensembl_marts, class = "data.frame")
      expect_identical(colnames(ensembl_marts), c("biomart", "version"))
    })
  },
  simplify = TRUE
)

test_that("Error when old URL is used", {
  expect_error(listMarts(host = "www.biomart.org"))
})
