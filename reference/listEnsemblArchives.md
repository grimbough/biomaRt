# Lists the available archived versions of Ensembl

Lists the available archived versions of Ensembl

## Usage

``` r
listEnsemblArchives(fetch = TRUE)
```

## Arguments

- fetch:

  If `TRUE` (default), fetches the list of Ensembl archives from the
  Ensembl website, and adds it to the `ensembl_versions` dataset. If
  `FALSE`, returns the existing `ensembl_versions` dataset.

  Returns a table containing the available archived versions of Ensembl,
  along with the dates they were created and the URL used to access
  them.

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
#> 30    Ensembl 104 May 2021 https://may2021.archive.ensembl.org     104
#> 31    Ensembl 103 Feb 2021 https://feb2021.archive.ensembl.org     103
#> 32    Ensembl 102 Nov 2020 https://nov2020.archive.ensembl.org     102
#> 33    Ensembl 101 Aug 2020 https://aug2020.archive.ensembl.org     101
#> 34    Ensembl 100 Apr 2020 https://apr2020.archive.ensembl.org     100
#> 35     Ensembl 99 Jan 2020 https://jan2020.archive.ensembl.org      99
#> 36     Ensembl 98 Sep 2019 https://sep2019.archive.ensembl.org      98
#> 37     Ensembl 97 Jul 2019 https://jul2019.archive.ensembl.org      97
#> 38     Ensembl 96 Apr 2019 https://apr2019.archive.ensembl.org      96
#> 39     Ensembl 95 Jan 2019 https://jan2019.archive.ensembl.org      95
#> 40     Ensembl 94 Oct 2018 https://oct2018.archive.ensembl.org      94
#> 41     Ensembl 93 Jul 2018 https://jul2018.archive.ensembl.org      93
#> 42     Ensembl 92 Apr 2018 https://apr2018.archive.ensembl.org      92
#> 43     Ensembl 91 Dec 2017 https://dec2017.archive.ensembl.org      91
#> 44     Ensembl 90 Aug 2017 https://aug2017.archive.ensembl.org      90
#> 45     Ensembl 89 May 2017 https://may2017.archive.ensembl.org      89
#> 46     Ensembl 88 Mar 2017 https://mar2017.archive.ensembl.org      88
#> 47     Ensembl 87 Dec 2016 https://dec2016.archive.ensembl.org      87
#> 48     Ensembl 86 Oct 2016 https://oct2016.archive.ensembl.org      86
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
#> 30                
#> 31                
#> 32                
#> 33                
#> 34                
#> 35                
#> 36                
#> 37                
#> 38                
#> 39                
#> 40                
#> 41                
#> 42                
#> 43                
#> 44                
#> 45                
#> 46                
#> 47                
#> 48                
#> 14                
#> 15                
#> 16                
#> 17                
```
