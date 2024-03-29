library(biomaRt) 

## example Mart object
ex_mart <- Mart(biomart = "ensembl", 
                dataset = "hsapiens_gene_ensembl",
                host = "www.ensembl.org",
                attributes = data.frame(
                  name = c("ensembl_gene_id", "gene_exon_intron", "hgnc_symbol"),
                  page = c("sequences", "sequences", "feature_page")
                ),
                filters = data.frame(
                  name = c("ensembl_gene_id", "hgnc_symbol", "gene_exon_intron"),
                  description = c("Ensembl ID", "HGNC Symbol", "Gene seq including exons and introns"),
                  type = c("id_list", "id_list", "text")
                )
)

test_that("useEnsembl() error handling is OK", { 
  
  expect_error(useEnsembl(), regexp = "You must provide the argument")
  
  expect_warning(useEnsembl(biomart = "genes", host = "https://ensembl.org"), 
                 regexp = "You cannot use the host")
  
})


with_mock_dir("all_500", {
  test_that("Ensembl mirror selection works", {
  
  expect_error(.chooseEnsemblMirror(mirror = NULL), 
               regexp = "Unable to query any Ensembl site") 
  
  })
},
simplify = FALSE)

with_mock_dir("www_OK", {
  test_that("Ensembl mirror selection works", {
    
    expect_message(.chooseEnsemblMirror(mirror = "useast"), regexp = "unresponsive") |>
      expect_equal("www")
    
  })
},
simplify = FALSE)

test_that("Ensembl URLs are constructed correctly", {
  
  ## no arguments ##
  expect_silent(.constructEnsemblURL()) |>
    expect_equal("https://www.ensembl.org")
  
  ## mirror ##
  expect_warning(.constructEnsemblURL(mirror = "INVALID_MIRROR"), 
                 regexp = "Invalid mirror\\. Select a mirror") %>%
    expect_equal("https://www.ensembl.org")
  
  expect_equal(.constructEnsemblURL(mirror = "useast"), 
               "https://useast.ensembl.org")
  
  ## GRCh ##
  expect_equal(.constructEnsemblURL(GRCh = 37), 
               "https://grch37.ensembl.org")
  
  expect_warning(.constructEnsemblURL(GRCh = 38), 
                 regex = "Only 37 can be specified for GRCh version") %>%
    expect_equal("https://www.ensembl.org")
  
  ## version ##
  expect_silent(.constructEnsemblURL(version = "100")) |>
               expect_equal("https://apr2020.archive.ensembl.org")
  
  expect_error(.constructEnsemblURL(version = "00"), 
               regexp = "Specified Ensembl version is not available")
  
  ## combinations ##
  expect_error(.constructEnsemblURL(version = "100", GRCh = 37), 
               regexp = "version or GRCh arguments cannot be used together")
  
  expect_warning(.constructEnsemblURL(mirror = "asia", version = 100), 
                 regexp = "version or GRCh arguments cannot be used together with the mirror argument") %>%
    expect_equal("https://apr2020.archive.ensembl.org")
  
  expect_warning(.constructEnsemblURL(mirror = "useast", GRCh = 37), 
                 regexp = "version or GRCh arguments cannot be used together with the mirror argument") %>%
    expect_equal("https://grch37.ensembl.org")
})

test_that("sequence correct code is used to get sequence based on ID type", {
  
  skip_if_not_installed('mockery')
  
  m <- mockery::mock(data.frame(hgnc_symbol = "STAT1", ensembl_gene_id = "ENSG00000115415"),
            data.frame(gene_exon_intron = "ACGTACGT", ensembl_gene_id = "ENSG00000115415"),
            data.frame(gene_exon_intron = "ACGTACGT", ensembl_gene_id = "ENSG00000115415"))
  
  mockery::stub(.getSequenceFromId, "getBM", m)

  expect_is(hgnc <- .getSequenceFromId(id = "STAT1", type = "hgnc_symbol", 
                               seqType = "gene_exon_intron", mart = ex_mart), 
            "data.frame")
  expect_is(ens_id <- .getSequenceFromId(id = "ENSG00000115415", type = "ensembl_gene_id", 
                                       seqType = "gene_exon_intron", mart = ex_mart), 
            "data.frame")
  expect_identical(hgnc$gene_exon_intron, ens_id$gene_exon_intron)

})

test_that("correct settings are applied for getSequence using coordinates.", {
  
  skip_if_not_installed('mockery')
  
  m <- mockery::mock(TRUE, cycle = TRUE)
  mockery::stub(.getSequenceFromCoords, "getBM", m)
  
  expect_true(.getSequenceFromCoords(chromosome = 1, start = 1, end = 2, type = "ensembl_gene_id", 
                                   seqType = "gene_exon_intron", mart = ex_mart))
  expect_true(.getSequenceFromCoords(chromosome = 1, start = 1, end = 2, upstream = 50,
                                     type = "ensembl_gene_id", seqType = "gene_exon_intron", 
                                     mart = ex_mart))
  expect_true(.getSequenceFromCoords(chromosome = 1, start = 1, end = 2, downstream = 50,
                                     type = "ensembl_gene_id", seqType = "gene_exon_intron", 
                                     mart = ex_mart))
  
  expect_length(mockery::mock_args(m)[[1]]$filters, 3L)
  expect_true("upstream_flank"   %in% names(mockery::mock_args(m)[[2]]$filters))
  expect_true("downstream_flank" %in% names(mockery::mock_args(m)[[3]]$filters))
  
  
})

test_that("getSequence deploys correct sub-function", {
  skip_if_not_installed('mockery')
  
  mockery::stub(getSequence, ".getSequenceFromCoords", TRUE)
  mockery::stub(getSequence, ".getSequenceFromId", FALSE)
  
  expect_true(getSequence(chromosome = 1, start = 1, end = 2, type = "ensembl_gene_id", 
                          seqType = "gene_exon_intron", mart = ex_mart))
  expect_false(getSequence(id = "STAT1", type = "hgnc_symbol", seqType = "gene_exon_intron", mart = ex_mart))
  
})

