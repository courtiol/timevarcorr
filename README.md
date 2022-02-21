
<!-- README.md is generated from README.Rmd. Please edit that file -->

# timevarcorr

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/fw)](https://CRAN.R-project.org/package=timevarcorr)
<!-- badges: end -->

This R package is at a very early stage of development. Many features
are missing and it may be buggy. Use at your own risk.

It aims at computing the correlation between 2 time-series. It performs
a non-parametric kernel smoothing (using a common bandwidth) of all
underlying componenents required for the computation of a correlation
coefficient (i.e. *x*, *y*, *x*<sup>2</sup>, *y*<sup>2</sup>,
*x* \* *y*). An automatic selection procedure for the bandwidth
parameter is implemented. Gaussian, box and Epanechnikov kernels can be
used. Both Pearson and Spearman correlation coefficients can be
estimated.

## Installation

You can install the development version of timevarcorr from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("courtiol/timevarcorr")
```

## Example

Simple example using base-R syntax:

``` r
library(timevarcorr)
d <- stockprice[1:500, ]
test <- with(d, tcor(x = SP500, y = FTSE100, t = DateID, kernel = "normal"))
#> [1] "h selected using LOO-CV = 60.9"
plot(test, type = "l", ylab = "Correlation", xlab = "Time", main = "SP500 vs FTSE100", las = 1)
```

<img src="man/figures/README-example-1.png" width="70%" />

Same example using tidyverse syntax:

``` r
library(tidyverse)
#> Warning in system("timedatectl", intern = TRUE): running command 'timedatectl'
#> had status 1
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
#> ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
#> ✓ tibble  3.1.4     ✓ dplyr   1.0.7
#> ✓ tidyr   1.2.0     ✓ stringr 1.4.0
#> ✓ readr   2.1.2     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()

stockprice |> 
  slice(1:500) |>
  summarise(tcor(x = SP500, y = FTSE100, t = DateID, kernel = "normal")) |>
  ggplot() +
    aes(x = t, y = r) +
    geom_line() +
    labs(title = "SP500 vs FTSE100", x = "Time", y = "Correlation") +
    theme_classic()
#> [1] "h selected using LOO-CV = 60.9"
```

<img src="man/figures/README-example2-1.png" width="70%" />

Alternative using `dplyr::mutate` showing gaps of observations in the
series:

``` r
stockprice |> 
  slice(1:500) |>
  mutate(r = tcor(x = SP500, y = FTSE100, t = DateID, kernel = "normal", keep.missing = TRUE)$r) |>
  ggplot() +
    aes(x = DateID, y = r) +
    geom_line() +
    labs(title = "SP500 vs FTSE100", x = "Time", y = "Correlation") +
    theme_classic()
#> [1] "h selected using LOO-CV = 60.9"
```

<img src="man/figures/README-example3-1.png" width="70%" />

You can do more. For exampe, you can use other kernels, fix the
bandwidth manually, or use the Spearman’s rather than the Pearson’s
correlation coefficient:

``` r
test2 <- with(d, tcor(x = SP500, y = FTSE100, t = DateID, cor.method = "spearman", kernel = "box", h = 10))
plot(test2, type = "l", ylab = "Correlation (Spearman)", xlab = "Time", las = 1)
```

<img src="man/figures/README-example4-1.png" width="70%" />

## Inspiration

The inspiration for this work was the reading of the following paper:

Choi, JE., Shin, D.W. Nonparametric estimation of time varying
correlation coefficient. J. Korean Stat. Soc. 50, 333–353 (2021).
<https://doi.org/10.1007/s42952-020-00073-6>
