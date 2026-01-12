# Display the analysis code from the 2009 Nature protocols paper

This function opens an editor displaying the analysis code of the Nature
Protocols 2009 paper

## Usage

``` r
NP2009code()
```

## Details

The [`edit()`](https://rdrr.io/r/utils/edit.html) function uses
`getOption("editor")` to select the editor. Use, for instance,
`options(editor="emacs")` to set another editor.

## See also

[`edit()`](https://rdrr.io/r/utils/edit.html)

## Author

Steffen Durinck, Wolfgang Huber

## Examples

``` r
if (FALSE) { # interactive()
NP2009code()
}
```
