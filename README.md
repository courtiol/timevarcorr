
<!-- README.md is generated from README.Rmd. Please edit that file -->

# timevarcorr

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/fw)](https://CRAN.R-project.org/package=timevarcorr)
<!-- badges: end -->

This R package is at an early stage of development. So, please only use
at your own risk.

It aims at computing the correlation between 2 time-series following the
method described in the following paper:

Choi, JE., Shin, D.W. Nonparametric estimation of time varying
correlation coefficient. J. Korean Stat. Soc. 50, 333–353 (2021).
<https://doi.org/10.1007/s42952-020-00073-6>

The chief idea is to perform a non-parametric kernel smoothing (using a
common bandwidth) of all underlying components required for the
computation of a correlation coefficient
(i.e. ![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x"),
![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y "y"),
![x^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%5E2 "x^2"),
![y^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y%5E2 "y^2"),
![x\*y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%2Ay "x*y")).

The automatic selection procedure for the bandwidth parameter proposed
in the paper is implemented in this package. The same goes for the
computation of confidence intervals.

We also implemented the possibility to use Epanechnikov, Gaussian, or
box kernels, as well as to estimate either the Pearson or the Spearman
correlation coefficient.

## Installation

You can install the development version of timevarcorr from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes") ## uncomment and run if you don't have this package installed
remotes::install_github("courtiol/timevarcorr")
```

That should suffice!

Note that this package relies so far on only one direct dependency –
[lpridge](https://github.com/cran/lpridge) – which itself depends on
nothing but a plain R install.

Nonetheless, in some of the examples below, I also rely on the
[tidyverse](https://github.com/tidyverse) ecosystem, so you would need
to install this as well to reproduce the content of this README:

``` r
install.packages("tidyverse")
```

## Examples

The main function of this package is called `tcor` and its documentation
is available here:

``` r
help(tcor, package = timevarcorr)
```

Very simple example using base-R syntax:

``` r
library(timevarcorr)
#> timevarcorr loaded; type ?tcor for help on this package.

d <- stockprice[1:500, ]
example1 <- with(d, tcor(x = SP500, y = FTSE100, t = DateID, kernel = "normal"))
#> 
#> You may set `nb.cores` to a number higher than 1 for faster computation.
#> [1] "h selected using LOO-CV = 60.9"
plot(example1, type = "l")
```

<img src="man/figures/README-example-1.png" width="70%" style="display: block; margin: auto;" />

Same example using tidyverse syntax (with confidence interval):

``` r
library(tidyverse)

d |> 
  summarise(tcor(x = SP500, y = FTSE100, t = DateID,
                 kernel = "normal", CI = TRUE)) |>
  ggplot() +
    aes(x = t, y = r, ymin = lwr, ymax = upr) +
    geom_ribbon(fill = "grey") +
    geom_line() +
    labs(title = "SP500 vs FTSE100", x = "Time", y = "Correlation") +
    theme_classic()
#> [1] "h selected using LOO-CV = 60.9"
```

<img src="man/figures/README-example2-1.png" width="70%" style="display: block; margin: auto;" />

Same example showing gaps of observations in the time series:

``` r
d |> 
  summarise(tcor(x = SP500, y = FTSE100, t = DateID,
                 kernel = "normal", CI = TRUE, keep.missing = TRUE)) |>
  ggplot() +
    aes(x = t, y = r, ymin = lwr, ymax = upr) +
    geom_ribbon(fill = "grey") +
    geom_line() +
    labs(title = "SP500 vs FTSE100", x = "Time", y = "Correlation") +
    theme_classic()
#> 
#> You may set `nb.cores` to a number higher than 1 for faster computation.
#> [1] "h selected using LOO-CV = 60.9"
```

<img src="man/figures/README-example3-1.png" width="70%" style="display: block; margin: auto;" />

You can do more. For example, you can use other kernels, fix the
bandwidth manually, or use the Spearman’s rather than the Pearson’s
correlation coefficient:

``` r
example2 <- with(d, tcor(x = SP500, y = FTSE100, t = DateID,
                 cor.method = "spearman", kernel = "box", h = 10))
plot(example2, type = "l")
```

<img src="man/figures/README-example4-1.png" width="70%" style="display: block; margin: auto;" />

You can also test the difference in correlation coefficients between two
time points:

``` r
example3 <- with(d, tcor(x = SP500, y = FTSE100, t = DateID, kernel = "normal", CI = TRUE))
#> 
#> You may set `nb.cores` to a number higher than 1 for faster computation.
#> [1] "h selected using LOO-CV = 60.9"
equality_test(example3, t1 = "2000-05-02", t2 = "2001-05-02")
#>           t1        r1         t2        r2   delta_r SE_delta_r   T_stat  df
#> 1 2000-05-02 0.4354486 2001-05-02 0.5721993 0.1367507  0.1224749 1.116561 910
#>           p
#> 1 0.2644768
```

Or you can test if specific time points (or all) differ from a reference
value:

``` r
ref_test(example3, t = c("2000-05-02", "2001-05-02"), r_ref = 0.5)
#>            t         r r_ref     delta_r SE_delta_r     T_stat  df         p
#> 1 2000-05-02 0.4354486   0.5 -0.06455140 0.10082726 -0.6402177 910 0.5221922
#> 2 2001-05-02 0.5721993   0.5  0.07219934 0.06952677  1.0384394 910 0.2993414
#>   p_adjustment
#> 1         none
#> 2         none
```

## Devel corner

README file compiled using `devtools::build_readme()`, with the
following setup:

``` r
devtools::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.1.3 (2022-03-10)
#>  os       Ubuntu 21.10
#>  system   x86_64, linux-gnu
#>  ui       X11
#>  language (EN)
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       Europe/Berlin
#>  date     2022-03-22
#>  pandoc   2.17.1.1 @ /usr/lib/rstudio/bin/quarto/bin/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package     * version    date (UTC) lib source
#>  assertthat    0.2.1      2019-03-21 [1] CRAN (R 4.1.3)
#>  backports     1.4.1      2021-12-13 [1] CRAN (R 4.1.3)
#>  brio          1.1.3      2021-11-30 [1] CRAN (R 4.1.3)
#>  broom         0.7.12     2022-01-28 [1] CRAN (R 4.1.3)
#>  cachem        1.0.6      2021-08-19 [1] CRAN (R 4.1.3)
#>  callr         3.7.0      2021-04-20 [1] CRAN (R 4.1.3)
#>  cellranger    1.1.0      2016-07-27 [1] CRAN (R 4.1.3)
#>  cli           3.2.0      2022-02-14 [1] CRAN (R 4.1.3)
#>  colorspace    2.0-3      2022-02-21 [1] CRAN (R 4.1.3)
#>  crayon        1.5.0      2022-02-14 [1] CRAN (R 4.1.3)
#>  DBI           1.1.2      2021-12-20 [1] CRAN (R 4.1.3)
#>  dbplyr        2.1.1      2021-04-06 [1] CRAN (R 4.1.3)
#>  desc          1.4.1      2022-03-06 [1] CRAN (R 4.1.3)
#>  devtools      2.4.3      2021-11-30 [1] CRAN (R 4.1.3)
#>  digest        0.6.29     2021-12-01 [1] CRAN (R 4.1.3)
#>  dplyr       * 1.0.8      2022-02-08 [1] CRAN (R 4.1.3)
#>  ellipsis      0.3.2      2021-04-29 [1] CRAN (R 4.1.3)
#>  evaluate      0.15       2022-02-18 [1] CRAN (R 4.1.3)
#>  fansi         1.0.2      2022-01-14 [1] CRAN (R 4.1.3)
#>  farver        2.1.0      2021-02-28 [1] CRAN (R 4.1.3)
#>  fastmap       1.1.0      2021-01-25 [1] CRAN (R 4.1.3)
#>  forcats     * 0.5.1      2021-01-27 [1] CRAN (R 4.1.3)
#>  fs            1.5.2      2021-12-08 [1] CRAN (R 4.1.3)
#>  generics      0.1.2      2022-01-31 [1] CRAN (R 4.1.3)
#>  ggplot2     * 3.3.5      2021-06-25 [1] CRAN (R 4.1.3)
#>  glue          1.6.2      2022-02-24 [1] CRAN (R 4.1.3)
#>  gtable        0.3.0      2019-03-25 [1] CRAN (R 4.1.3)
#>  haven         2.4.3      2021-08-04 [1] CRAN (R 4.1.3)
#>  highr         0.9        2021-04-16 [1] CRAN (R 4.1.3)
#>  hms           1.1.1      2021-09-26 [1] CRAN (R 4.1.3)
#>  htmltools     0.5.2      2021-08-25 [1] CRAN (R 4.1.3)
#>  httr          1.4.2      2020-07-20 [1] CRAN (R 4.1.3)
#>  jsonlite      1.8.0      2022-02-22 [1] CRAN (R 4.1.3)
#>  knitr         1.37       2021-12-16 [1] CRAN (R 4.1.3)
#>  labeling      0.4.2      2020-10-20 [1] CRAN (R 4.1.3)
#>  lifecycle     1.0.1      2021-09-24 [1] CRAN (R 4.1.3)
#>  lubridate     1.8.0      2021-10-07 [1] CRAN (R 4.1.3)
#>  magrittr      2.0.2      2022-01-26 [1] CRAN (R 4.1.3)
#>  memoise       2.0.1      2021-11-26 [1] CRAN (R 4.1.3)
#>  modelr        0.1.8      2020-05-19 [1] CRAN (R 4.1.3)
#>  munsell       0.5.0      2018-06-12 [1] CRAN (R 4.1.3)
#>  pillar        1.7.0      2022-02-01 [1] CRAN (R 4.1.3)
#>  pkgbuild      1.3.1      2021-12-20 [1] CRAN (R 4.1.3)
#>  pkgconfig     2.0.3      2019-09-22 [1] CRAN (R 4.1.3)
#>  pkgload       1.2.4      2021-11-30 [1] CRAN (R 4.1.3)
#>  prettyunits   1.1.1      2020-01-24 [1] CRAN (R 4.1.3)
#>  processx      3.5.2      2021-04-30 [1] CRAN (R 4.1.3)
#>  ps            1.6.0      2021-02-28 [1] CRAN (R 4.1.3)
#>  purrr       * 0.3.4      2020-04-17 [1] CRAN (R 4.1.3)
#>  R6            2.5.1      2021-08-19 [1] CRAN (R 4.1.3)
#>  Rcpp          1.0.8.3    2022-03-17 [1] CRAN (R 4.1.3)
#>  readr       * 2.1.2      2022-01-30 [1] CRAN (R 4.1.3)
#>  readxl        1.3.1      2019-03-13 [1] CRAN (R 4.1.3)
#>  remotes       2.4.2      2021-11-30 [1] CRAN (R 4.1.3)
#>  reprex        2.0.1      2021-08-05 [1] CRAN (R 4.1.3)
#>  rlang         1.0.2      2022-03-04 [1] CRAN (R 4.1.3)
#>  rmarkdown     2.13       2022-03-10 [1] CRAN (R 4.1.3)
#>  rprojroot     2.0.2      2020-11-15 [1] CRAN (R 4.1.3)
#>  rstudioapi    0.13       2020-11-12 [1] CRAN (R 4.1.3)
#>  rvest         1.0.2      2021-10-16 [1] CRAN (R 4.1.3)
#>  scales        1.1.1      2020-05-11 [1] CRAN (R 4.1.3)
#>  sessioninfo   1.2.2      2021-12-06 [1] CRAN (R 4.1.3)
#>  stringi       1.7.6      2021-11-29 [1] CRAN (R 4.1.3)
#>  stringr     * 1.4.0      2019-02-10 [1] CRAN (R 4.1.3)
#>  testthat      3.1.2      2022-01-20 [1] CRAN (R 4.1.3)
#>  tibble      * 3.1.6      2021-11-07 [1] CRAN (R 4.1.3)
#>  tidyr       * 1.2.0      2022-02-01 [1] CRAN (R 4.1.3)
#>  tidyselect    1.1.2      2022-02-21 [1] CRAN (R 4.1.3)
#>  tidyverse   * 1.3.1      2021-04-15 [1] CRAN (R 4.1.3)
#>  timevarcorr * 0.0.0.9005 2022-03-22 [1] local
#>  tzdb          0.2.0      2021-10-27 [1] CRAN (R 4.1.3)
#>  usethis       2.1.5      2021-12-09 [1] CRAN (R 4.1.3)
#>  utf8          1.2.2      2021-07-24 [1] CRAN (R 4.1.3)
#>  vctrs         0.3.8      2021-04-29 [1] CRAN (R 4.1.3)
#>  withr         2.5.0      2022-03-03 [1] CRAN (R 4.1.3)
#>  xfun          0.30       2022-03-02 [1] CRAN (R 4.1.3)
#>  xml2          1.3.3      2021-11-30 [1] CRAN (R 4.1.3)
#>  yaml          2.3.5      2022-02-21 [1] CRAN (R 4.1.3)
#> 
#>  [1] /home/courtiol/R/x86_64-pc-linux-gnu-library/4.1
#>  [2] /usr/local/lib/R/site-library
#>  [3] /usr/lib/R/site-library
#>  [4] /usr/lib/R/library
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
