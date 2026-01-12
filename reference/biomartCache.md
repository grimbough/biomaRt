# biomaRt result caching

biomaRt makes use of a results cache to speedup execution of queries
that have been run before. These functions provide details on the status
of this cache, and allow it to be deleted.

## Usage

``` r
biomartCacheClear()

biomartCacheInfo()
```

## Value

These functions do not return anything and are called for their side
effects. `biomartCacheInfo()` prints the location of the cache, along
with the number of files and their total size on disk.
`biomartCacheClear()` will delete the current contents of the cache.

## Author

Mike Smith
