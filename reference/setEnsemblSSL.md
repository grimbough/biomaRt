# Save system specific SSL settings for contacting Ensembl

On some systems specific SSL settings have to be applied to allow https
connections to the Ensembl servers. This function allows these to be
saved in the biomaRt cache, so they will be retrieved each time they are
needed. biomaRt will try to determine them automatically, but this
function can be used to set them manually if required.

## Usage

``` r
setEnsemblSSL(settings)
```

## Arguments

- settings:

  A named list. Each entry should be a valid curl option, as found in
  [`curl::curl_options()`](https://jeroen.r-universe.dev/curl/reference/curl_options.html).

## Author

Mike Smith

## Examples

``` r
if (FALSE) { # \dontrun{
ssl_settings <- list(
  "ssl_cipher_list" = "DEFAULT@SECLEVEL=1",
  "ssl_verifypeer"  = FALSE
)
setEnsemblSSL(ssl_settings)
} # }
```
