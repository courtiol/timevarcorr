#' Stockprice
#'
#' A dataset containing the stockmarket returns between 2000-04-03 and 2017-12-05.
#' This dataset is very close to the one used by Choi & Shin 2021, although not
#' strictly identical. It has been produced by the Oxford-Man Institute of Quantitative Finance.
#'
#' @format A data frame with 4618 rows and 7 variables:
#' \describe{
#'   \item{DateID}{a vector of `Date`.}
#'   \item{SP500}{a numeric vector of the stockmarket return for the S&P 500 Index.}
#'   \item{FTSE100}{a numeric vector of the stockmarket return for the FTSE 100.}
#'   \item{Nikkei}{a numeric vector of the stockmarket return for the Nikkei 225.}
#'   \item{DAX}{a numeric vector of the stockmarket return for the German stock index.}
#'   \item{NASDAQ}{a numeric vector of the stockmarket return for the Nasdaq Stock Market.}
#'   \item{Event}{a character string of particular events that have impacted the stockmarket, as in Choi & Shin 2021.}
#' }
#'
#' @source
#' - General webpage:
#' \url{https://realized.oxford-man.ox.ac.uk/data/download}
#'
#' - Direct download link:
#' \url{https://realized.oxford-man.ox.ac.uk/images/oxfordmanrealizedvolatilityindices-0.2-final.zip}
#'
#' @seealso
#' [`datasets::EuStockMarkets`] for a similar dataset, albeit formated differently.
#'
#' @references
#' Heber, Gerd, Asger Lunde, Neil Shephard and Kevin Sheppard (2009)
#' "Oxford-Man Institute's realized library", Oxford-Man Institute, University of Oxford.
#'
#' Choi, JE., Shin, D.W. Nonparametric estimation of time varying correlation coefficient.
#' J. Korean Stat. Soc. 50, 333–353 (2021). https://doi.org/10.1007/s42952-020-00073-6
#'
"stockprice"
