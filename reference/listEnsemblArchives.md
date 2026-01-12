# Lists the available archived versions of Ensembl

Returns a table containing the available archived versions of Ensembl,
along with the dates they were created and the URL used to access them.

## Usage

``` r
listEnsemblArchives()
```

## Author

Mike Smith

## Examples

``` r
listEnsemblArchives()
#> Error in req_perform(html_request): Failed to perform HTTP request.
#> Caused by error in `curl::curl_fetch_memory()`:
#> ! Timeout was reached [www.ensembl.org]:
#> Operation timed out after 10002 milliseconds with 0 bytes received
```
