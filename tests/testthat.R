library(testthat)
library(mockery)
library(biomaRt)

biomartCacheClear()

test_check("biomaRt", encoding = "UTF-8")