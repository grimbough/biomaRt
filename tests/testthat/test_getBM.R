library(biomaRt)
context('Testing methods getBM() function')

expect_error(getBM(), "You must provide a valid Mart object")