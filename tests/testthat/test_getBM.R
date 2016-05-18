library(biomaRt)
context('Testing methods getBM() function')

expect_error(getBM(), "Error in martCheck\\(mart\\)")