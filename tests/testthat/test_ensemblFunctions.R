library(biomaRt)

test_that("useEnsembl() error handling is OK", { 
  
  expect_error(useEnsembl(), regexp = "You must provide the argument")
  
})

test_that("Ensembl URLs are constructed correctly", {
  
  ## no arguments ##
  expect_equal(.constructEnsemblURL(),
               "https://www.ensembl.org")
  
  ## mirror ##
  expect_warning(.constructEnsemblURL(mirror = "INVALID_MIRROR"), 
                 regexp = "Invalid mirror\\. Select a mirror") %>%
    expect_equal("https://www.ensembl.org")
  
  expect_equal(.constructEnsemblURL(mirror = "uswest"), 
               "https://uswest.ensembl.org")
  
  ## GRCh ##
  expect_equal(.constructEnsemblURL(GRCh = 37), 
               "https://grch37.ensembl.org")
  
  expect_warning(.constructEnsemblURL(GRCh = 38), 
                 regex = "Only 37 can be specified for GRCh version") %>%
    expect_equal("https://www.ensembl.org")
  
  ## version ##
  expect_equal(.constructEnsemblURL(version = "100"), 
               "https://apr2020.archive.ensembl.org")
  
  expect_error(.constructEnsemblURL(version = "00"), 
               regexp = "Specified Ensembl version is not available")
  
  ## combinations ##
  expect_error(.constructEnsemblURL(version = "100", GRCh = 37), 
               regexp = "version or GRCh arguments cannot be used together")
  
  expect_warning(.constructEnsemblURL(mirror = "asia", version = 100), 
                 regexp = "version or GRCh arguments cannot be used together with the mirror argument") %>%
    expect_equal("https://apr2020.archive.ensembl.org")
  
  expect_warning(.constructEnsemblURL(mirror = "uswest", GRCh = 37), 
                 regexp = "version or GRCh arguments cannot be used together with the mirror argument") %>%
    expect_equal("https://grch37.ensembl.org")
})
