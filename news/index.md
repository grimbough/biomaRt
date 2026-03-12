# Changelog

## CHANGES IN biomaRt VERSION 2.68.0

### BREAKING CHANGES

- The `port` argument, if unspecified and if no scheme is provided in
  the `host` argument, now defaults to `443` (the standard port for
  `https`), instead of `80` (the standard port for `http`). This fixes
  issue [\#155](https://github.com/Huber-group-EMBL/biomaRt/issues/155)
  reported by Lori Shepherd.

### MINOR IMPROVEMENTS

- All requests sent by biomaRt now include a custom user agent linking
  to the source GitHub repository, and indicating the package version.

## CHANGES IN biomaRt VERSION 2.66.0

### BREAKING CHANGES

- the `https` argument in
  [`useEnsembl()`](https://huber-group-embl.github.io/biomaRt/reference/useEnsembl.md)
  has been removed, as announced in version 2.50.0. Connections to
  Ensembl mirrors now always use `https://`, as enforced by a redirect
  from the server.

### USER VISIBLE CHANGES

- This package documentation is now available as a pkgdown website at
  <https://huber-group-embl.github.io/biomaRt/>.

- The
  [`exportFASTA()`](https://huber-group-embl.github.io/biomaRt/reference/exportFASTA.md)
  function is now roughly 5x faster.

### BUG FIXES

- [`getGene()`](https://huber-group-embl.github.io/biomaRt/reference/getGene.md)
  no longer errors with “This function only works when used with the
  ensembl BioMart.”, that was always returned, even when using the
  ensembl BioMart. This bug was introduced in version 2.27.1
  (Bioconductor 3.3). Thanks to [Tobias
  Hoch](https://github.com/toobiwankenobi) for the report.

- This package now explicitly depends on R \>= 4.1.0. This was
  effectively already the case since version 2.58, in line with
  Bioconductor policy, but not explicitly documented in `DESCRIPTION`.

- [`getBM()`](https://huber-group-embl.github.io/biomaRt/reference/getBM.md)
  no longer silently convert alleles “T” and “F” to logical `TRUE` and
  `FALSE`. This rarely happened when a single allele, coincidentally
  named “T” or “F”, was returned in a column. Thanks to [Sina
  Rüeger](https://github.com/sinarueeger) for the report.

- [`getBM()`](https://huber-group-embl.github.io/biomaRt/reference/getBM.md)
  now internally handles redirects and automatically resubmits the query
  to the new location. This fixes issues when using versioned forms such
  as <https://eg59-plants.ensembl.org/>, which redirects to
  <https://may2024-plants.ensembl.org/index.html>. Thanks to Edoardo
  Bertolini for the report.

### INTERNAL CHANGES

- Coding style throughout the package has been harmonized using the air
  tool. Contributors using RStudio, Positron or VS Code should have
  their code styled automatically on save.

- A duplicated definition of the `.getEnsemblSSL()` internal function
  has been removed.

- Static analysis via the lintr package is now performed on each push
  and PR. It should mostly be invisible to users but might result in
  slightly increased performance in some cases.

- This package now uses roxygen2 to generate documentation and
  `NAMESPACE`.

- A `R CMD check` `NOTE` about a missing import has been resolved.

- This package continuous integration setup now errors on `R CMD check`
  `WARNING`s. This is checked on each push and PR.

- The digest and rappdirs dependencies have been removed in favour of
  functions provided by the base R tools package.

- Examples have been reviewed and fixed where necessary.

- testthat edition 3 (instead of previously edition 2) is now used for
  unit tests.

## CHANGES IN biomaRt VERSION 2.64.0

### BUG FIXES

- Version 1.10 of {httr2} changed how URLs are parsed, and this broke
  some biomaRt functionality. This has been patched. (Backported to
  biomaRt 2.62.1)

## CHANGES IN biomaRt VERSION 2.62.0

### USER VISIBLE CHANGES

- Several deprecated functions are now defunct.

### BUG FIXES

- Results returned from BioMart queries will be read using Latin-1
  encounding if the default fails. Reported in
  <https://support.bioconductor.org/p/9158844/> (Backported to 2.60.1)

- Fixed issue when only one dataset was listed in a Mart instance,
  causing data.frame dimensions to be dropped. This broke connectivity
  to <https://parasite.wormbase.org>. (Backported to 2.60.1)

## CHANGES IN biomaRt VERSION 2.60.0

### USER VISIBLE CHANGES

- listEnsemblGenomes() and useEnsemblGenomes() now have a host argument,
  allowing you to select an Ensembl Genomes archive site. (Thanks to
  Hervé Pagès [@hpages](https://github.com/hpages) for the suggestion:
  <https://github.com/grimbough/biomaRt/issues/93>)

- The ‘curl’ argument to getBM() has been deprecated as it is no longer
  applicable and doesn’t do anything.

### INTERNAL CHANGES

- Removed dependency on XML package and switched all functionality to
  xml2

- Swiched from httr to httr2 package for submitting queries to BioMart
  servers.

## CHANGES IN biomaRt VERSION 2.58.0

### USER VISIBLE CHANGES

- getSequence() will now provide a more informative error message if
  requesting a flanking sequence and not provided with an upstream or
  downstream range.

- Remove references to the uswest mirror, which has now been retired
  (<https://www.ensembl.info/2023/01/13/retirement-of-ensembl-us-west-aws-mirror/>)

## CHANGES IN biomaRt VERSION 2.56.0

### BUG FIXES

- Fix problem when multiple cache entries with the same ID could be
  created. (Thanks to Hervé Pagès & Henrik Bengtsson for independent
  reports of this issue.)

- bmRequest() will now respect the setting in options(“timeout”)

## CHANGES IN biomaRt VERSION 2.52.0

### BUG FIXES

- Stop reporting message about the use of https when using useEnsembl()
  with a ‘version’ argument.

- Use virtualSchemaName provided by a Mart, rather than simply
  “default”. This caused issues with the Ensembl Plants Mart.

## CHANGES IN biomaRt VERSION 2.50.0

### MINOR CHANGES

- useMart() and listMarts() will warn users if using http to access
  Ensembl. https will be enforced by Ensembl from late 2021.

### BUG FIXES

- Address issue where checking the list of Ensembl Archives would stop
  all queries from working if the main www.ensembl.org site was
  unavailable.

- Fix bug introduced in getSequence() where asking for flanking
  sequences resulted in an invalid query.

- The argument ‘host’ is no longer ignored in useEnsembl() (Thanks to
  forum user “A” - <https://support.bioconductor.org/p/9139019/>)

## CHANGES IN biomaRt VERSION 2.48.0

### NEW FEATURES

- getSequence() now allows the cache to be turned off via the ‘useCache’
  argument.

- Automatic detection of SSL issues with Ensembl, and appropriate
  settings applied to httr functions used by biomaRt.

### BUG FIXES

- Addressed issue with getSequence() and ID types that are not available
  on the ‘sequences’ page. This could result in truncated sequences
  being returned from a query.

- getBM() would fail if it found a cache entry, but the file was
  corrupted. Invalid entries are now detected and deleted if
  encountered.

## CHANGES IN biomaRt VERSION 2.46.0

### BUG FIXES

- getLDS() now detects if trying to use datasets from different Marts
  and reports this to the user.

## CHANGES IN biomaRt VERSION 2.42.0

### NEW FEATURES

- The results of queries will now be cached, and if repeated queries are
  detected the results are loaded from disk.

### MINOR CHANGES

- Ensembl users will be redirected to their closest mirror unless the
  host argument is explicitly provided. In this case the defined value
  will be enforced.

- Unused argument ‘ssl.verifypeer’ removed from listMarts() and
  useMarts().

- RCurl removed from package dependecies.

### BUG FIXES

- Improvements made to selecting the correct port when using http vs
  https

- Results that contain unescaped new line characters are now returned
  successfully.

## CHANGES IN biomaRt VERSION 2.36.0

### BUG FIXES

- Patched problem returning the list of available datasets, if the
  description of one or more datasets included an apostrophe (introduced
  with new primate species in Ensembl).

- Caught scenario where ensemblRedirect=FALSE was still being ignored.

- Changed query submission when redirection is detected to cope with
  apparently new behaviour of the Ensembl mirrors.

### MINOR CHANGES

- Increase query timeout limit to 5 minutes.

## CHANGES IN biomaRt VERSION 2.34.0

### NEW FEATURES

- Added the listEnsemblArchives() function. This returns a table of the
  available Ensembl archives, and replaces the archive = TRUE argument
  to several functions, which was no longer working.

### BUG FIXES

- The Ensembl BioMart server doesn’t always respond well if queries with
  more than 500 filter values are submitted. If a query that exceed this
  is detect biomaRt will now submit the query in batches and concatonate
  the result when completed.

### MINOR CHANGES

- You can now provide a host with ‘<http://>’ at the start, or a
  trailing ‘/’ (typically copy/pasted from a browser) and useMarts() etc
  will cope.

## CHANGES IN biomaRt VERSION 2.32.0

### BUG FIXES

- Fixed bug when columns were not returned in the order requested, which
  resulted in the wrong column names being added to the result.

## CHANGES IN biomaRt VERSION 2.30.0

### SIGNIFICANT USER-LEVEL CHANGES

- Updated vignette to use BiocStyle and execute most code chunks.
