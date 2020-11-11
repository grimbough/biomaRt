
test_that("Both SSL settings are applied if needed", { 
  
  ## clear any httr configuration before running these tests
  httr::reset_config()

  m <- mock(stop("sslv3 alert handshake failure"),
            stop("unable to get local issuer certificate"))
  stub(.check_ensembl_ssl, ".test_ensembl", m)
  
  expect_silent(.check_ensembl_ssl())
  expect_length(getOption("httr_config")$options, 2L)
  expect_identical(getOption("httr_config")$options$ssl_cipher_list, "DEFAULT@SECLEVEL=1")
  expect_false(getOption("httr_config")$options$ssl_verifypeer)
  
})

test_that("Only one SSL setting is applied", { 

  ## clear any httr configuration before running these tests
  httr::reset_config()
  
  m <- mock(stop("sslv3 alert handshake failure"),
            "",
            stop("unable to get local issuer certificate"),
            "")
  stub(.check_ensembl_ssl, ".test_ensembl", m)
  
  expect_silent(.check_ensembl_ssl())
  expect_length(getOption("httr_config")$options, 1L)
  expect_identical(getOption("httr_config")$options$ssl_cipher_list, "DEFAULT@SECLEVEL=1")
  
  ## clear again
  httr::reset_config()
  
  expect_silent(.check_ensembl_ssl())
  expect_length(getOption("httr_config")$options, 1L)
  expect_false(getOption("httr_config")$options$ssl_verifypeer)
  
})
  



