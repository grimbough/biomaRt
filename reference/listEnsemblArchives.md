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
#>              name     date                                 url version
#> 1  Ensembl GRCh37 Feb 2014          https://grch37.ensembl.org  GRCh37
#> 2     Ensembl 116 Jun 2026 https://jun2026.archive.ensembl.org     116
#> 3     Ensembl 115 Sep 2025 https://sep2025.archive.ensembl.org     115
#> 4     Ensembl 114 May 2025 https://may2025.archive.ensembl.org     114
#> 5     Ensembl 113 Oct 2024 https://oct2024.archive.ensembl.org     113
#> 6     Ensembl 112 May 2024 https://may2024.archive.ensembl.org     112
#> 7     Ensembl 111 Jan 2024 https://jan2024.archive.ensembl.org     111
#> 8     Ensembl 110 Jul 2023 https://jul2023.archive.ensembl.org     110
#> 9     Ensembl 109 Feb 2023 https://feb2023.archive.ensembl.org     109
#> 10    Ensembl 108 Oct 2022 https://oct2022.archive.ensembl.org     108
#> 11    Ensembl 107 Jul 2022 https://jul2022.archive.ensembl.org     107
#> 12    Ensembl 106 Apr 2022 https://apr2022.archive.ensembl.org     106
#> 13    Ensembl 105 Dec 2021 https://dec2021.archive.ensembl.org     105
#> 14     Ensembl 80 May 2015 https://may2015.archive.ensembl.org      80
#> 15     Ensembl 77 Oct 2014 https://oct2014.archive.ensembl.org      77
#> 16     Ensembl 75 Feb 2014 https://feb2014.archive.ensembl.org      75
#> 17     Ensembl 54 May 2009 https://may2009.archive.ensembl.org      54
#>    current_release
#> 1                 
#> 2                *
#> 3                 
#> 4                 
#> 5                 
#> 6                 
#> 7                 
#> 8                 
#> 9                 
#> 10                
#> 11                
#> 12                
#> 13                
#> 14                
#> 15                
#> 16                
#> 17                
```
