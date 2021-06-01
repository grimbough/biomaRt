

test_that("Archive list fails if all mirrors are giving error page", { 
  
  skip_if_not_installed('mockery')
  
  ## mock the output from httr::GET to mimic the ensembl error page
  mock_get <- mockery::mock(
    httr:::response(url = "www.ensembl.org", 
                    status_code = 200, 
                    content = charToRaw("The Ensembl web service you requested is temporarily unavailable")),
    cycle = TRUE
  )
  mockery::stub(.getArchiveList, "GET", mock_get)

  expect_error(.getArchiveList(), regexp = "Unable to contact any Ensembl mirror")
})

test_that("Archive list fails if all mirrors don't give 200 responses", { 
  
  skip_if_not_installed('mockery')
  
  ## mock the output from httr::GET to mimic the ensembl error page
  mock_get <- mockery::mock(
    httr:::response(url = "www.ensembl.org", 
                    status_code = 400, 
                    content = charToRaw("This is an error")),
    cycle = TRUE
  )
  mockery::stub(.getArchiveList, "GET", mock_get)
  
  expect_error(.getArchiveList(), regexp = "Unable to contact any Ensembl mirror")
})
