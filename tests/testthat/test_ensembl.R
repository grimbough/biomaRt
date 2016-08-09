library(biomaRt)
context('Testing ensembl specific functions')

## invalid mirror
expect_warning(mirror <- .useEnsemblMirror("NOT_A_MIRROR"), "Invalid mirror")
expect_identical(mirror, "www.ensembl.org")
## valid mirrors
expect_equal(.useEnsemblMirror("uswest"), "uswest.ensembl.org")
expect_equal(.useEnsemblMirror("useast"), "useast.ensembl.org")
expect_equal(.useEnsemblMirror("asia"), "asia.ensembl.org")

## requesting specific versions
expect_equal(.useEnsemblVersion(version = 84), "e84.ensembl.org")
expect_equal(.useEnsemblVersion(GRCh = 37), "grch37.ensembl.org")
## can only specifc GRCh37, anthing else warns and uses default
expect_warning(host <- .useEnsemblVersion(GRCh = 38), "Only 37 can be specified for GRCh version")
expect_identical(host, "www.ensembl.org")
