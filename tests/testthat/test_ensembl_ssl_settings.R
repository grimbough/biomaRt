cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)
bfc <- BiocFileCache::BiocFileCache(biomartCacheInfo(), ask = FALSE)

test_that("We can manually set SSL settings", {
  ## empty list
  settings <- list()
  expect_true(setEnsemblSSL(settings = settings))
  expect_identical(settings, .getEnsemblSSL())

  ## one element
  settings <- list("foo" = "A")
  expect_true(setEnsemblSSL(settings = settings))
  expect_identical(settings, .getEnsemblSSL())

  ## add element
  settings <- list("baa" = "B")
  expect_true(setEnsemblSSL(settings = settings))
  expect_identical(list("foo" = "A", "baa" = "B"), .getEnsemblSSL())

  ## update existing element
  settings <- list("baa" = "C")
  expect_true(setEnsemblSSL(settings = settings))
  expect_identical(list("foo" = "A", "baa" = "C"), .getEnsemblSSL())

  ## remove these entries from the cache so they aren't used in later tests
  biomartCacheClear()
})
