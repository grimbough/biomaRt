# Using a BioMart other than Ensembl

## Introduction

In recent years a wealth of biological data has become available in
public data repositories. Easy access to these valuable data resources
and firm integration with data analysis is needed for comprehensive
bioinformatics data analysis. The
*[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* package,
provides an interface to a growing collection of databases implementing
the [BioMart software
suite](https://www.ensembl.org/info/data/biomart/index.html). The
package enables retrieval of large amounts of data in a uniform way
without the need to know the underlying database schemas or write
complex SQL queries. Examples of BioMart databases are Ensembl, Uniprot
and HapMap. These major databases give
*[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* users direct
access to a diverse set of data and enable a wide range of powerful
online queries from R.

## Using a BioMart other than Ensembl

There are a small number of non-Ensembl databases that offer a BioMart
interface to their data. The
*[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* package can
be used to access these in a very similar fashion to Ensembl. The
majority of *[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)*
functions will work in the same manner, but the construction of the
initial Mart object requires slightly more setup. In this section we
demonstrate the setting requires to query [Wormbase
ParaSite](https://parasite.wormbase.org/index.html) and
[Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html). First we need
to load *[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)*.

``` r
library(biomaRt)
```

### Wormbase

To demonstrate the use of the
*[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* package with
non-Ensembl databases the next query is performed using the Wormbase
ParaSite BioMart. In this example, we use the
[`listMarts()`](https://huber-group-embl.github.io/biomaRt/reference/listMarts.md)
function to find the name of the available marts, given the URL of
Wormbase. We use this to connect to Wormbase BioMart using the
[`useMart()`](https://huber-group-embl.github.io/biomaRt/reference/useMart.md)
function.[¹](#fn1)

``` r
listMarts(host = "parasite.wormbase.org")
```

    ## Error in `req_perform()`:
    ## ! Failed to perform HTTP request.
    ## Caused by error in `curl::curl_fetch_memory()`:
    ## ! Unsupported protocol [parasite.wormbase.org]:
    ## Received HTTP/0.9 when not allowed

``` r
wormbase <- useMart(
  biomart = "parasite_mart",
  host = "https://parasite.wormbase.org",
  port = 443
)
```

    ## Error in `req_perform()`:
    ## ! HTTP 500 Internal Server Error.

We can then use functions described earlier in this vignette to find and
select the gene dataset, and print the first 6 available attributes and
filters. Then we use a list of gene names as filter and retrieve
associated transcript IDs and the transcript biotype.

``` r
listDatasets(wormbase)
```

    ## Error:
    ## ! object 'wormbase' not found

``` r
wormbase <- useDataset(mart = wormbase, dataset = "wbps_gene")
```

    ## Error:
    ## ! object 'wormbase' not found

``` r
head(listFilters(wormbase))
```

    ## Error:
    ## ! object 'wormbase' not found

``` r
head(listAttributes(wormbase))
```

    ## Error:
    ## ! object 'wormbase' not found

``` r
getBM(
  attributes = c(
    "external_gene_id",
    "wbps_transcript_id",
    "transcript_biotype"
  ),
  filters = "gene_name",
  values = c("unc-26", "his-33"),
  mart = wormbase
)
```

    ## Error:
    ## ! object 'wormbase' not found

### Phytozome

#### Version 12

The Phytozome 12 BioMart was
[retired](https://jgi.doe.gov/more-intuitive-phytozome-interface/) in
August 2021 and can not longer be accessed.

#### Version 13

Version 13 of Phytozome can be found at
<https://phytozome-next.jgi.doe.gov/> and if you wish to query that
version the URL used to create the Mart object must reflect that.

``` r
phytozome_v13 <- useMart(
  biomart = "phytozome_mart",
  dataset = "phytozome",
  host = "https://phytozome-next.jgi.doe.gov"
)
```

Once this is set up the usual
*[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* functions
can be used to interrogate the database options and run queries.

``` r
getBM(
  attributes = c("organism_name", "gene_name1"),
  filters = "gene_name_filter",
  values = "82092",
  mart = phytozome_v13
)
```

    ##          organism_name gene_name1
    ## 1 Smoellendorffii_v1.0      82092

## Session Info

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8    
    ##  [5] LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C             
    ##  [9] LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] biomaRt_2.67.6   BiocStyle_2.38.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] KEGGREST_1.50.0      xfun_0.57            bslib_0.10.0         httr2_1.2.2         
    ##  [5] Biobase_2.70.0       vctrs_0.7.2          tools_4.5.3          generics_0.1.4      
    ##  [9] stats4_4.5.3         curl_7.0.0           tibble_3.3.1         AnnotationDbi_1.72.0
    ## [13] RSQLite_2.4.6        blob_1.3.0           pkgconfig_2.0.3      dbplyr_2.5.2        
    ## [17] desc_1.4.3           S4Vectors_0.48.0     lifecycle_1.0.5      compiler_4.5.3      
    ## [21] stringr_1.6.0        textshaping_1.0.5    Biostrings_2.78.0    progress_1.2.3      
    ## [25] Seqinfo_1.0.0        htmltools_0.5.9      sass_0.4.10          yaml_2.3.12         
    ## [29] pillar_1.11.1        pkgdown_2.2.0        crayon_1.5.3         jquerylib_0.1.4     
    ## [33] cachem_1.1.0         tidyselect_1.2.1     digest_0.6.39        stringi_1.8.7       
    ## [37] purrr_1.2.1          dplyr_1.2.0          bookdown_0.46        fastmap_1.2.0       
    ## [41] cli_3.6.5            magrittr_2.0.4       withr_3.0.2          prettyunits_1.2.0   
    ## [45] filelock_1.0.3       rappdirs_0.3.4       bit64_4.6.0-1        rmarkdown_2.31      
    ## [49] XVector_0.50.0       httr_1.4.8           bit_4.6.0            ragg_1.5.2          
    ## [53] png_0.1-9            hms_1.1.4            memoise_2.0.1        evaluate_1.0.5      
    ## [57] knitr_1.51           IRanges_2.44.0       BiocFileCache_3.0.0  rlang_1.1.7         
    ## [61] glue_1.8.0           DBI_1.3.0            BiocManager_1.30.27  xml2_1.5.2          
    ## [65] BiocGenerics_0.56.0  jsonlite_2.0.0       R6_2.6.1             systemfonts_1.3.2   
    ## [69] fs_2.0.1

``` r
warnings()
```

------------------------------------------------------------------------

1.  Note that we use the `https` address and must provide the port as
    `443`. Queries to WormBase will fail without these options.
