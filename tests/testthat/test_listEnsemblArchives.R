
with_mock_dir("all_err_page", {
  test_that("Archive list fails if all mirrors are giving error page", {
    
    expect_error(.getArchiveList(), regexp = "Unable to contact any Ensembl mirror")
  })
}, 
simplify = TRUE)

with_mock_dir("all_500", {
  test_that("Archive list fails if all mirrors don't give 200 responses", {
  
    expect_error(.getArchiveList(), regexp = "Unable to contact any Ensembl mirror")
  })
}, 
simplify = TRUE)

