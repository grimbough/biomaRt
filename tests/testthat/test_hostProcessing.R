library(biomaRt)

## adding http if needed
host <- 'www.myurl.org'
expect_equal(object = .cleanHostURL(host = host),
             expected = "http://www.myurl.org")
 
## stripping trailing slash
host <- 'http://www.myurl.org/'
expect_equal(object = .cleanHostURL(host = host),
             expected = "http://www.myurl.org")

## leave https already there
host <- 'https://www.myurl.org'
expect_equal(object = .cleanHostURL(host = host),
             expected = "https://www.myurl.org")
