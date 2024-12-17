
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

gopher_pcoa[1:5, 1:5]
#>             [,1]         [,2]        [,3]          [,4]         [,5]
#> [1,]  0.06099355 -0.097959065 -0.06764845 -0.0070942187 -0.008393500
#> [2,]  0.05965239 -0.094230934 -0.06432898 -0.0061189343 -0.007090959
#> [3,] -0.10934852  0.008355788 -0.01432736 -0.0050424583 -0.008221191
#> [4,] -0.08247512  0.003757215 -0.00715867 -0.0007547401  0.001563042
#> [5,] -0.11000241  0.008495579 -0.01444203 -0.0051544261 -0.008313776
lice_pcoa[1:5, 1:5]
#>             [,1]        [,2]         [,3]          [,4]          [,5]
#> [1,] -0.01582840 -0.09886551  0.069847344 -0.0027606488 -0.1152976632
#> [2,] -0.01348795 -0.07749434  0.044488641 -0.0027442380 -0.0627717149
#> [3,] -0.08906150  0.04633471 -0.001014638 -0.0008995698  0.0009282499
#> [4,] -0.12890429  0.08997653  0.006351793  0.0003770117 -0.0038739941
#> [5,] -0.08406230  0.04226708 -0.001770707 -0.0011051837  0.0009239243
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
#> 1     1        2 0.0009312199 0.004004004
#> 2     1        8 0.0011611181 0.094094094
#> 3     2        1 0.0010501927 0.014014014
#> 4     2        9 0.0006711532 0.113113113
#> 5     3        3 0.0017178115 0.001001001
#> 6     4        7 0.0010412160 0.017017017
#> 7     5        5 0.0016337474 0.001001001
#> 8     6        4 0.0019256509 0.002002002
#> 9     7        6 0.0015769140 0.006006006
#> 10    8       16 0.0012866890 0.014014014
#> 11    9       14 0.0007419528 0.181181181
#> 12   10       17 0.0013933007 0.014014014
#> 13   11       15 0.0015419676 0.001001001
#> 14   12       10 0.0007215688 0.310310310
#> 15   13       11 0.0004458883 0.587587588
#> 16   14       13 0.0006382121 0.371371371
#> 17   15       12 0.0008493063 0.071071071
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
#> 1     1        2 0.0009312199 0.006006006
#> 2     1        8 0.0011611181 0.098098098
#> 3     2        1 0.0010501927 0.020020020
#> 4     2        9 0.0006711532 0.109109109
#> 5     3        3 0.0017178115 0.001001001
#> 6     4        7 0.0010412160 0.021021021
#> 7     5        5 0.0016337474 0.001001001
#> 8     6        4 0.0019256509 0.008008008
#> 9     7        6 0.0015769140 0.009009009
#> 10    8       16 0.0012866890 0.018018018
#> 11    9       14 0.0007419528 0.199199199
#> 12   10       17 0.0013933007 0.012012012
#> 13   11       15 0.0015419676 0.009009009
#> 14   12       10 0.0007215688 0.298298298
#> 15   13       11 0.0004458883 0.592592593
#> 16   14       13 0.0006382121 0.393393393
#> 17   15       12 0.0008493063 0.070070070
```
