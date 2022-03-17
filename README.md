
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
computation of a correlation coefficient (i.e. *x*, *y*,
*x*<sup>2</sup>, *y*<sup>2</sup>, *x* \* *y*).

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

You can also test the difference in correlation coefficents between two
time points:

``` r
example3 <- with(d, tcor(x = SP500, y = FTSE100, t = DateID, kernel = "normal", CI = TRUE))
#> 
#> You may set `nb.cores` to a number higher than 1 for faster computation.
#> [1] "h selected using LOO-CV = 60.9"
equality_test(example3, t1 = "2000-05-02", t2 = "2001-05-02")
#>     delta_r SE_delta_r   T_stat df_student pv_student
#> 1 0.1367507  0.1224749 1.116561        910  0.2644768
```

## Devel corner

README file compiled using `devtools::build_readme()`, with the
following setup:

``` r
devtools::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.1.2 (2021-11-01)
#>  os       Fedora Linux 37 (Container Image Prerelease)
#>  system   x86_64, linux-gnu
#>  ui       X11
#>  language (EN)
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       Europe/Berlin
#>  date     2022-03-17
#>  pandoc   2.14.0.3 @ /usr/libexec/rstudio/bin/pandoc/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package     * version    date (UTC) lib source
#>  assertthat    0.2.1      2019-03-21 [2] CRAN (R 4.1.0)
#>  backports     1.4.1      2021-12-13 [2] CRAN (R 4.1.2)
#>  brio          1.1.3      2021-11-30 [2] CRAN (R 4.1.2)
#>  broom         0.7.12     2022-01-28 [2] CRAN (R 4.1.2)
#>  cachem        1.0.6      2021-08-19 [2] CRAN (R 4.1.0)
#>  callr         3.7.0      2021-04-20 [2] CRAN (R 4.1.0)
#>  cellranger    1.1.0      2016-07-27 [2] CRAN (R 4.1.0)
#>  cli           3.1.1      2022-01-20 [2] CRAN (R 4.1.2)
#>  colorspace    2.0-2      2021-06-24 [2] CRAN (R 4.1.0)
#>  CoprManager   0.3.9      2022-01-04 [4] local
#>  crayon        1.4.2      2021-10-29 [2] CRAN (R 4.1.1)
#>  DBI           1.1.1      2021-01-15 [2] CRAN (R 4.1.0)
#>  dbplyr        2.1.1      2021-04-06 [2] CRAN (R 4.1.0)
#>  desc          1.4.0      2021-09-28 [2] CRAN (R 4.1.1)
#>  devtools      2.4.2      2021-06-07 [2] CRAN (R 4.1.0)
#>  digest        0.6.28     2021-09-23 [2] CRAN (R 4.1.1)
#>  dplyr       * 1.0.7      2021-06-18 [2] CRAN (R 4.1.0)
#>  ellipsis      0.3.2      2021-04-29 [2] CRAN (R 4.1.0)
#>  evaluate      0.14       2019-05-28 [2] CRAN (R 4.1.0)
#>  fansi         0.5.0      2021-05-25 [2] CRAN (R 4.1.0)
#>  farver        2.1.0      2021-02-28 [2] CRAN (R 4.1.0)
#>  fastmap       1.1.0      2021-01-25 [2] CRAN (R 4.1.0)
#>  forcats     * 0.5.1      2021-01-27 [2] CRAN (R 4.1.0)
#>  fs            1.5.2      2021-12-08 [2] CRAN (R 4.1.2)
#>  generics      0.1.2      2022-01-31 [2] CRAN (R 4.1.2)
#>  ggplot2     * 3.3.5      2021-06-25 [2] CRAN (R 4.1.0)
#>  glue          1.6.1      2022-01-22 [2] CRAN (R 4.1.2)
#>  gtable        0.3.0      2019-03-25 [2] CRAN (R 4.1.0)
#>  haven         2.4.3      2021-08-04 [2] CRAN (R 4.1.0)
#>  highr         0.9        2021-04-16 [2] CRAN (R 4.1.0)
#>  hms           1.1.1      2021-09-26 [2] CRAN (R 4.1.1)
#>  htmltools     0.5.2      2021-08-25 [2] CRAN (R 4.1.1)
#>  httr          1.4.2      2020-07-20 [2] CRAN (R 4.1.0)
#>  jsonlite      1.7.2      2020-12-09 [2] CRAN (R 4.1.0)
#>  knitr         1.36       2021-09-29 [2] CRAN (R 4.1.1)
#>  labeling      0.4.2      2020-10-20 [2] CRAN (R 4.1.0)
#>  lifecycle     1.0.1      2021-09-24 [2] CRAN (R 4.1.1)
#>  lubridate     1.8.0      2021-10-07 [2] CRAN (R 4.1.1)
#>  magrittr      2.0.2      2022-01-26 [2] CRAN (R 4.1.2)
#>  memoise       2.0.1      2021-11-26 [2] CRAN (R 4.1.2)
#>  modelr        0.1.8      2020-05-19 [2] CRAN (R 4.1.0)
#>  munsell       0.5.0      2018-06-12 [2] CRAN (R 4.1.0)
#>  pillar        1.7.0      2022-02-01 [2] CRAN (R 4.1.2)
#>  pkgbuild      1.3.0      2021-12-09 [2] CRAN (R 4.1.2)
#>  pkgconfig     2.0.3      2019-09-22 [2] CRAN (R 4.1.0)
#>  pkgload       1.2.4      2021-11-30 [2] CRAN (R 4.1.2)
#>  prettyunits   1.1.1      2020-01-24 [2] CRAN (R 4.1.0)
#>  processx      3.5.2      2021-04-30 [2] CRAN (R 4.1.0)
#>  ps            1.6.0      2021-02-28 [2] CRAN (R 4.1.0)
#>  purrr       * 0.3.4      2020-04-17 [2] CRAN (R 4.1.0)
#>  R6            2.5.1      2021-08-19 [2] CRAN (R 4.1.0)
#>  Rcpp          1.0.7      2021-07-07 [2] CRAN (R 4.1.0)
#>  readr       * 2.1.2      2022-01-30 [2] CRAN (R 4.1.2)
#>  readxl        1.3.1      2019-03-13 [2] CRAN (R 4.1.0)
#>  remotes       2.4.2      2021-11-30 [2] CRAN (R 4.1.2)
#>  reprex        2.0.1      2021-08-05 [2] CRAN (R 4.1.0)
#>  rlang         1.0.0      2022-01-26 [2] CRAN (R 4.1.2)
#>  rmarkdown     2.11       2021-09-14 [2] CRAN (R 4.1.1)
#>  rprojroot     2.0.2      2020-11-15 [2] CRAN (R 4.1.0)
#>  rstudioapi    0.13       2020-11-12 [2] CRAN (R 4.1.0)
#>  rvest         1.0.2      2021-10-16 [2] CRAN (R 4.1.1)
#>  scales        1.1.1      2020-05-11 [2] CRAN (R 4.1.0)
#>  sessioninfo   1.2.2      2021-12-06 [2] CRAN (R 4.1.2)
#>  stringi       1.7.4      2021-08-25 [2] CRAN (R 4.1.1)
#>  stringr     * 1.4.0      2019-02-10 [2] CRAN (R 4.1.0)
#>  testthat      3.1.2      2022-01-20 [2] CRAN (R 4.1.2)
#>  tibble      * 3.1.4      2021-08-25 [2] CRAN (R 4.1.1)
#>  tidyr       * 1.2.0      2022-02-01 [2] CRAN (R 4.1.2)
#>  tidyselect    1.1.1      2021-04-30 [2] CRAN (R 4.1.0)
#>  tidyverse   * 1.3.1      2021-04-15 [2] CRAN (R 4.1.0)
#>  timevarcorr * 0.0.0.9003 2022-03-17 [1] local
#>  tzdb          0.2.0      2021-10-27 [2] CRAN (R 4.1.1)
#>  usethis       2.1.5      2021-12-09 [2] CRAN (R 4.1.2)
#>  utf8          1.2.2      2021-07-24 [2] CRAN (R 4.1.0)
#>  vctrs         0.3.8      2021-04-29 [2] CRAN (R 4.1.0)
#>  withr         2.4.3      2021-11-30 [2] CRAN (R 4.1.2)
#>  xfun          0.27       2021-10-18 [2] CRAN (R 4.1.1)
#>  xml2          1.3.3      2021-11-30 [2] CRAN (R 4.1.2)
#>  yaml          2.2.2      2022-01-25 [2] CRAN (R 4.1.2)
#> 
#>  [1] /home/alex/R/x86_64-redhat-linux-gnu-library/4.1
#>  [2] /usr/local/lib/R/library
#>  [3] /usr/lib64/R/library
#>  [4] /usr/share/R/library
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
