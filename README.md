
<!-- README.md is generated from README.Rmd. Please edit that file -->

# parafit

<!-- badges: start -->

![GitHub R package
version](https://img.shields.io/github/r-package/v/shaunson26/parafit)
<!-- badges: end -->

The goal of parafit is to …

## Installation

You can install the development version of parafit from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Shaunson26/parafit")
```

Functions of the `parafit` package have the suffix `pf_` e.g.

``` r
pf_pcoa
pf_parafit
```

    #> ℹ Loading parafit

## Example

Create Principal Coordinates of phylogentic distances.

``` r
# gopher.D and lice.D are included in parafit and are square distance matrices
gopher_pcoa <- pf_pcoa(gopher.D)
lice_pcoa <- pf_pcoa(lice.D)

gopher_pcoa[1:5, 1:3]
#>             [,1]         [,2]        [,3]
#> [1,]  0.06099355 -0.097959065 -0.06764845
#> [2,]  0.05965239 -0.094230934 -0.06432898
#> [3,] -0.10934852  0.008355788 -0.01432736
#> [4,] -0.08247512  0.003757215 -0.00715867
#> [5,] -0.11000241  0.008495579 -0.01444203
lice_pcoa[1:5, 1:3]
#>             [,1]        [,2]         [,3]
#> [1,] -0.01582840 -0.09886551  0.069847344
#> [2,] -0.01348795 -0.07749434  0.044488641
#> [3,] -0.08906150  0.04633471 -0.001014638
#> [4,] -0.12890429  0.08997653  0.006351793
#> [5,] -0.08406230  0.04226708 -0.001770707
```

Run parafit obtaining the global test statistics and testing individual
associations

``` r
res <-
  pf_parafit(host_pcoa = t(gopher_pcoa),
             parasite_pcoa = lice_pcoa,
             associations = gopher.lice.links,
             permutations = 999,
             test_links = TRUE)
#> Running parafit ...
#> There are 17 links to be tested
#> Done

res
#> Parafit
#> 
#> ParaFitGlobal = 0.01389872, P value = 0.001 (999)
#> 
#> Individual tests of links
#> 
#>    host parasite       stat_1         p_1
#> 1     1        2 0.0009312199 0.007007007
#> 2     1        8 0.0011611181 0.106106106
#> 3     2        1 0.0010501927 0.020020020
#> 4     2        9 0.0006711532 0.107107107
#> 5     3        3 0.0017178115 0.001001001
#> 6     4        7 0.0010412160 0.019019019
#> 7     5        5 0.0016337474 0.001001001
#> 8     6        4 0.0019256509 0.008008008
#> 9     7        6 0.0015769140 0.011011011
#> 10    8       16 0.0012866890 0.018018018
#> 11    9       14 0.0007419528 0.195195195
#> 12   10       17 0.0013933007 0.021021021
#> 13   11       15 0.0015419676 0.006006006
#> 14   12       10 0.0007215688 0.313313313
#> 15   13       11 0.0004458883 0.591591592
#> 16   14       13 0.0006382121 0.384384384
#> 17   15       12 0.0008493063 0.054054054
```

For larger datasets, testing individual associations can be run in
parallel

``` r
res <-
  pf_parafit(host_pcoa = t(gopher_pcoa),
             parasite_pcoa = lice_pcoa,
             associations = gopher.lice.links,
             permutations = 999,
             test_links = TRUE,
             parallel = TRUE,
             cores = 2)
#> Running parafit ...
#> There are 17 links to be tested ... in parallel
#> Done

res
#> Parafit
#> 
#> ParaFitGlobal = 0.01389872, P value = 0.001 (999)
#> 
#> Individual tests of links
#> 
#>    host parasite       stat_1         p_1
#> 1     1        2 0.0009312199 0.009009009
#> 2     1        8 0.0011611181 0.098098098
#> 3     2        1 0.0010501927 0.022022022
#> 4     2        9 0.0006711532 0.107107107
#> 5     3        3 0.0017178115 0.001001001
#> 6     4        7 0.0010412160 0.015015015
#> 7     5        5 0.0016337474 0.001001001
#> 8     6        4 0.0019256509 0.003003003
#> 9     7        6 0.0015769140 0.005005005
#> 10    8       16 0.0012866890 0.016016016
#> 11    9       14 0.0007419528 0.177177177
#> 12   10       17 0.0013933007 0.023023023
#> 13   11       15 0.0015419676 0.005005005
#> 14   12       10 0.0007215688 0.293293293
#> 15   13       11 0.0004458883 0.600600601
#> 16   14       13 0.0006382121 0.385385385
#> 17   15       12 0.0008493063 0.066066066
```
