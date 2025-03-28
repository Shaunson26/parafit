
<!-- README.md is generated from README.Rmd. Please edit that file -->

# parafit

<!-- badges: start -->

<figure>
<img src="https://img.shields.io/github/r-package/v/shaunson26/parafit"
alt="GitHub R package version" />
<figcaption aria-hidden="true">GitHub R package version</figcaption>
</figure>

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

## Example

Create Principal Coordinates of phylogenetic distances. Coordinates link
the vectors element of the output.

``` r
# gopher.D and lice.D are included in parafit and are square distance matrices
gopher_pcoa <- pf_pcoa(gopher.D)
lice_pcoa <- pf_pcoa(lice.D)

str(gopher_pcoa)
#> List of 3
#>  $ pco       :List of 5
#>   ..$ correction: chr [1:2] "none" "1"
#>   ..$ note      : chr "There were no negative eigenvalues. No correction was applied"
#>   ..$ values    :'data.frame':   14 obs. of  5 variables:
#>   .. ..$ Eigenvalues   : num [1:14] 0.0731 0.0384 0.0329 0.0236 0.016 ...
#>   .. ..$ Relative_eig  : num [1:14] 0.3043 0.1597 0.1369 0.0984 0.0666 ...
#>   .. ..$ Broken_stick  : num [1:14] 0.2323 0.1608 0.1251 0.1013 0.0834 ...
#>   .. ..$ Cumul_eig     : num [1:14] 0.304 0.464 0.601 0.699 0.766 ...
#>   .. ..$ Cumul_br_stick: num [1:14] 0.232 0.393 0.518 0.619 0.703 ...
#>   ..$ vectors   : num [1:15, 1:14] 0.061 0.0597 -0.1093 -0.0825 -0.11 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:15] "T.talpoides" "T.bottae" "O.underwoodi" "O.hispidus" ...
#>   .. .. ..$ : chr [1:14] "Axis.1" "Axis.2" "Axis.3" "Axis.4" ...
#>   ..$ trace     : num 0.24
#>   ..- attr(*, "class")= chr "pcoa"
#>  $ vectors   : num [1:15, 1:14] 0.061 0.0597 -0.1093 -0.0825 -0.11 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:15] "T.talpoides" "T.bottae" "O.underwoodi" "O.hispidus" ...
#>   .. ..$ : chr [1:14] "Axis.1" "Axis.2" "Axis.3" "Axis.4" ...
#>  $ sum_eig_sq: num 0.00918

gopher_pcoa$vectors[1:5, 1:3]
#>                   Axis.1       Axis.2      Axis.3
#> T.talpoides   0.06099355  0.097959065 -0.06764845
#> T.bottae      0.05965239  0.094230934 -0.06432898
#> O.underwoodi -0.10934852 -0.008355788 -0.01432736
#> O.hispidus   -0.08247512 -0.003757215 -0.00715867
#> O.cavator    -0.11000241 -0.008495579 -0.01444203
lice_pcoa$vectors[1:5, 1:3]
#>                   Axis.1      Axis.2       Axis.3
#> T.minor      -0.01582840 -0.09886551  0.069847344
#> T.barbarae   -0.01348795 -0.07749434  0.044488641
#> G.setzeri    -0.08906150  0.04633471 -0.001014638
#> G.cherriei   -0.12890429  0.08997653  0.006351793
#> G.panamensis -0.08406230  0.04226708 -0.001770707
```

Run parafit obtaining the global test statistics and testing individual
associations

``` r
res <-
  pf_parafit(host_pcoa = gopher_pcoa$vectors,
             parasite_pcoa = lice_pcoa$vectors,
             associations = gopher.lice.links,
             permutations = 999,
             test_links = TRUE)
#> Running parafit ...
#> There are 17 links to be tested
#> Testing links ...
#> 1
#> 17
#> Done

res
#> Parafit
#> 
#> ParaFitGlobal = 0.01389872, P value = 0.001 (999)
#> 
#> Individual tests of links
#> 
#>    host parasite       stat_1         p_1
#> 1     1        2 0.0009312199 0.008008008
#> 2     1        8 0.0011611181 0.084084084
#> 3     2        1 0.0010501927 0.020020020
#> 4     2        9 0.0006711532 0.092092092
#> 5     3        3 0.0017178115 0.001001001
#> 6     4        7 0.0010412160 0.011011011
#> 7     5        5 0.0016337474 0.001001001
#> 8     6        4 0.0019256509 0.003003003
#> 9     7        6 0.0015769140 0.003003003
#> 10    8       16 0.0012866890 0.022022022
#> 11    9       14 0.0007419528 0.185185185
#> 12   10       17 0.0013933007 0.013013013
#> 13   11       15 0.0015419676 0.003003003
#> 14   12       10 0.0007215688 0.293293293
#> 15   13       11 0.0004458883 0.609609610
#> 16   14       13 0.0006382121 0.357357357
#> 17   15       12 0.0008493063 0.049049049
```

For larger datasets, testing individual associations can be run in
parallel

``` r
res <-
  pf_parafit(host_pcoa = gopher_pcoa$vectors,
             parasite_pcoa = lice_pcoa$vectors,
             associations = gopher.lice.links,
             permutations = 999,
             test_links = TRUE,
             parallel = TRUE,
             cores = 2)
#> Running parafit ...
#> There are 17 links to be tested ... in parallel
#> Testing links ...
#> 1
#> 17
#> Done

res
#> Parafit
#> 
#> ParaFitGlobal = 0.01389872, P value = 0.001 (999)
#> 
#> Individual tests of links
#> 
#>    host parasite       stat_1         p_1
#> 1     1        2 0.0009312199 0.005005005
#> 2     1        8 0.0011611181 0.091091091
#> 3     2        1 0.0010501927 0.008008008
#> 4     2        9 0.0006711532 0.109109109
#> 5     3        3 0.0017178115 0.001001001
#> 6     4        7 0.0010412160 0.007007007
#> 7     5        5 0.0016337474 0.001001001
#> 8     6        4 0.0019256509 0.004004004
#> 9     7        6 0.0015769140 0.005005005
#> 10    8       16 0.0012866890 0.018018018
#> 11    9       14 0.0007419528 0.185185185
#> 12   10       17 0.0013933007 0.019019019
#> 13   11       15 0.0015419676 0.006006006
#> 14   12       10 0.0007215688 0.287287287
#> 15   13       11 0.0004458883 0.603603604
#> 16   14       13 0.0006382121 0.388388388
#> 17   15       12 0.0008493063 0.067067067
```
