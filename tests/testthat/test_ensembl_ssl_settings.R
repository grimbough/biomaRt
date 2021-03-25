cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)
bfc <- BiocFileCache::BiocFileCache(biomartCacheInfo(), ask = FALSE)

test_that("Both SSL settings are applied if needed", { 
    
    skip_if_not_installed('mockery')
    
    ## clear any httr configuration before running these tests
    httr::reset_config()
    
    m <- mockery::mock(stop("sslv3 alert handshake failure"),
              stop("unable to get local issuer certificate"), 
              "")
    mockery::stub(.checkEnsemblSSL, ".test_ensembl", m)
    
    expect_silent(httr_config <- .checkEnsemblSSL())
    expect_is(httr_config, "list")
    expect_length(httr_config, 2L)
    expect_identical(httr_config[[1]]$options$ssl_cipher_list, "DEFAULT@SECLEVEL=1")
    expect_false(httr_config[[2]]$options$ssl_verifypeer)
    
})

test_that("Only one SSL setting is applied", { 

    skip_if_not_installed('mockery')

    m <- mockery::mock(stop("sslv3 alert handshake failure"),
              "",
              stop("unable to get local issuer certificate"),
              "")
    mockery::stub(.checkEnsemblSSL, ".test_ensembl", m)
    
    expect_silent(httr_config <- .checkEnsemblSSL())
    expect_length(httr_config, 1L)
    expect_identical(httr_config[[1]]$options$ssl_cipher_list, "DEFAULT@SECLEVEL=1")
    
    expect_silent(httr_config <- .checkEnsemblSSL())
    expect_length(httr_config, 1L)
    expect_false(httr_config[[1]]$options$ssl_verifypeer)
    
})

test_that("If we hit an unknown error", { 
    
    skip_if_not_installed('mockery')
    
    ## clear any httr configuration before running these tests
    httr::reset_config()
    
    m <- mockery::mock(stop("This is an unexpected SSL error"), 
              "")
    mockery::stub(.checkEnsemblSSL, ".test_ensembl", m)
    
    expect_message(.checkEnsemblSSL(), "Possible SSL connectivity problems detected")
    expect_length(getOption("httr_config")$options, 0L)
    
})


test_that("SSL settings are stored in the cache", {
    
    skip_if_not_installed('mockery')
    
    m <- mockery::mock(stop("sslv3 alert handshake failure"),
                       stop("unable to get local issuer certificate"), 
                       "")
    mockery::stub(.checkEnsemblSSL, ".test_ensembl", m)
    
    expect_silent(httr_config <- .getEnsemblSSL())
    expect_is(httr_config, "list")
    expect_true(.checkInCache(bfc, "ensembl-ssl-settings"))
})
