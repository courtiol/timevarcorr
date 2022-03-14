#' Compute time varying correlation coefficients
#'
#' The function `tcor` implements the nonparametric estimation of the time
#' varying correlation coefficient proposed by Choi & Shin, 2021. The general
#' idea is to compute a (Pearson) correlation coefficient (`r(x,y) = (mean(xy) -
#' mean(x)*mean(y)) / (sqrt(mean(x^2)-mean(x)^2) * sqrt(mean(y^2)-mean(y)^2))`),
#' but instead of using the means required for such a computation, each
#' component (i.e. `x`, `y`, `x^2`, `y^2`, `x*y`) is smoothed and the smoothed
#' terms are considered in place the original means. The intensity of the
#' smoothing depends on a unique parameter: the bandwidth (`h`). If `h = Inf`,
#' the method produces the original (i.e. time-invariant) correlation value. The
#' smaller the parameter `h`, the more variation in time is being captured. The
#' parameter `h` can be provided by the user; otherwise it is automatically
#' estimated (see **Details**).
#'
#' - **Smoothing**: the smoothing of each component is performed by kernel
#' regression. The default is to use the Epanechnikov kernel following Choi &
#' Shin 2021, but other kernels have also been implemented and can thus
#' alternatively be used (see [`kern_smooth`] for details). The normal kernell
#' seems to sometimes lead to very small bandwidth being selected, but the
#' default can lead to numerical issues. We thus recommend always comparing the
#' results from different kernel methods.
#'
#' - **Bandwidth selection**: when the value used to define the bandwidth (`h`)
#' is set to `NULL` (the default), it is first estimated by leave-one-out cross
#' validation. If cross validation error is minimal for the maximal value of `h`
#' considered (`8*sqrt(N)`), rather than taking this as the optimal `h` value,
#' the bandwidth becomes estimated using the so-called elbow criterion. This
#' latter method identifies the value `h` after which the cross validation error
#' decreasing very little. The procedure is detailed in section 2.1 in Choi &
#' Shin, 2021.
#'
#' - **Parallel computation**: if `h` is not provided, an automatic bandwidth
#' selection occurs (see above). For large datasets, this step can be
#' computationally demanding. The current implementation thus relies on
#' [`parallel::mclapply`] on thus is only effective for Linux and MacOS. Relying
#' on parallel processing also implies that you call `options("mc.cores" = XX)`
#' beforehand, replacing `XX` by the relevant number of CPU cores you want to use
#' (see **Examples**).
#'
#' @inheritParams kern_smooth
#' @param x a numeric vector.
#' @param y a numeric vector of to be correlated with `x`.
#' @param cor.method a character string indicating which correlation coefficient
#' is to be computed ("pearson", the default; or "spearman").
#' @param keep.missing a logical specifying if time steps associated with missing
#' information should be kept in the output (default = `FALSE` to facilitate plotting).
#' @param verbose a logical specifying if information should be displayed to
#' monitor the progress of the cross validation (default = `FALSE`).
#'
#' @return A dataframe of time points (`t`) and corresponding correlation values (`r`).
#' Some metadata are also attached to the dataframe (as attributes):
#' - `h` the bandwidth parameter.
#' - `CV_error` the minimal Cross Validation error when `h` selected by CV.
#' - `h_selection` the method used to select `h`.
#' - `h_select_duration` the computing time spent to select the bandwidth
#' parameter.
#'
#' @references
#' Choi, JE., Shin, D.W. Nonparametric estimation of time varying correlation coefficient.
#' J. Korean Stat. Soc. 50, 333–353 (2021). https://doi.org/10.1007/s42952-020-00073-6
#'
#' @seealso [`kern_smooth`]
#'
#' @importFrom parallel mclapply
#' @export
#'
#' @examples
#'
#' ## Effect of the bandwidth
#'
#' res_h50  <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = 50))
#' res_h100  <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = 100))
#' res_h200 <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = 200))
#' plot(res_h50, type = "l", ylab = "Cor", xlab = "Time", las = 1, col = "grey")
#' points(res_h100, type = "l", col = "blue")
#' points(res_h200, type = "l", col = "red")
#' legend("topright", fill = c("grey", "blue", "red"),
#'        legend = c("50", "100", "200"), bty = "n", title = "Bandwidth (h)")
#'
#'
#' ## Effect of the correlation method
#'
#' res_pearson <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = 150))
#' res_spearman <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = 150,
#'                                      cor.method = "spearman"))
#' plot(res_pearson, type = "l", ylab = "Cor", xlab = "Time", las = 1)
#' points(res_spearman, type = "l", col = "blue")
#' legend("topright", fill = c("black", "blue"),
#'        legend = c("pearson", "spearman"), bty = "n", title = "cor.method")
#'
#'
#' ## Infinite bandwidth and fixed correlation correspondance
#' ## nb: does not work with default kernel
#'
#' res_pearson_hInf  <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = Inf,
#'                                            kernel = "normal"))
#' res_spearman_hInf <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = Inf,
#'                                            kernel = "normal", cor.method = "spearman"))
#' r <- cor(stockprice$SP500, stockprice$FTSE100, use = "pairwise.complete.obs")
#' rho <- cor(stockprice$SP500, stockprice$FTSE100, method = "spearman", use = "pairwise.complete.obs")
#' plot(res_pearson_hInf, type = "l", ylim = c(0, 1), ylab = "Cor", xlab = "Time", las = 1, lwd = 4)
#' points(res_spearman_hInf, type = "l", lwd = 4, col = "blue")
#' abline(h = r, col = "red", lty = 2, lwd = 2)
#' abline(h = rho, col = "green", lty = 2, lwd = 2)
#'
#'
#' \dontrun{
#' ## Automatic selection of the bandwidth using parallel processing and comparison
#' ## of the 3 alternative kernels on full dataset
#' # nb: takes a few minutes to run
#'
#' options("mc.cores" = 2L)
#' res_hauto_epanech <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID,
#'                                           verbose = TRUE))
#' res_hauto_box <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID,
#'                                        kernel = "box", verbose = TRUE))
#' res_hauto_norm <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID,
#'                                         kernel = "norm", verbose = TRUE))
#' plot(res_hauto_epanech, type = "l", col = "red",
#'      ylab = "Cor", xlab = "Time", las = 1, ylim = c(0, 1))
#' points(res_hauto_box, type = "l", col = "blue")
#' points(res_hauto_norm, type = "l", col = "orange")
#' legend("topright", fill = c("red", "blue", "orange"),
#'        legend = c("epanechnikov", "box", "normal"), bty = "n",
#'        title = "Kernel")
#'
#'
#' ## Cross validation error according to each kernel
#'
#' attr(res_hauto_epanech, "CV_error")
#' attr(res_hauto_box, "CV_error")
#' attr(res_hauto_norm, "CV_error")
#'
#' ## Selected bandwidth according to each kernel
#'
#' attr(res_hauto_epanech, "h")
#' attr(res_hauto_box, "h")
#' attr(res_hauto_norm, "h")
#'
#'
#' ## Comparison of the 3 alternative kernels under same bandwidth
#'
#' res_epanech <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID,
#'                                      h = attr(res_hauto_epanech, "h")))
#' res_box <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, kernel = "box",
#'                                  h = attr(res_hauto_epanech, "h")))
#' res_norm <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, kernel = "norm",
#'                                   h = attr(res_hauto_epanech, "h")))
#'
#' plot(res_epanech, type = "l", col = "red", ylab = "Cor", xlab = "Time",
#'      las = 1, ylim = c(0, 1))
#' points(res_box, type = "l", col = "grey")
#' points(res_norm, type = "l", col = "orange")
#' legend("topright", fill = c("red", "grey", "orange"),
#'        legend = c("epanechnikov", "box", "normal"), bty = "n",
#'        title = "Kernel")
#'
#'
#' ## Not all kernels work well in all situations
#' ## nb: EuStockMarkets is a time-series object provided with R
#'
#' EuStock_epanech <- tcor(EuStockMarkets[, "DAX"], EuStockMarkets[, "SMI"])
#' plot(EuStock_epanech, type = "l", las = 1) ## problem of estimation for first and last years
#'
#' EuStock_norm <- tcor(EuStockMarkets[, "DAX"], EuStockMarkets[, "SMI"], kernel = "normal")
#' plot(EuStock_norm, type = "l", las = 1) ## normal kernel seems to work great with these data
#'
#' }
#'
tcor <- function(x, y, t = seq_along(x), h = NULL, cor.method = c("pearson", "spearman"),
                 kernel = c("epanechnikov", "box", "normal"), param_smoother = list(),
                 keep.missing = FALSE,
                 verbose = FALSE) {

  ## stripping out missing data
  missing <- is.na(x) | is.na(y) | is.na(t)
  x_ori <- x[!missing] ## ori = original, i.e. not smoothed
  y_ori <- y[!missing]
  t_ori <- t[!missing]

  ## if the bandwidth is not given, we must estimate it
  CV_error <- NA # default value for export if not computed
  if (is.null(h)) {

    ## check nb.cores input
    if (is.null(options("mc.cores")[[1]])) {
      nb.cores <- 1
    } else {
      nb.cores <- options("mc.cores")[[1]]
    }

    if (Sys.info()[['sysname']] != "Windows" && parallel::detectCores() < nb.cores) {
      stop(paste("\nYour computer does not allow such a large value for the argument `nb.cores`. The maximum value you may consider is", parallel::detectCores(), ".\n"))
    }

    if (Sys.info()[['sysname']] != "Windows" && parallel::detectCores() > 1 && nb.cores == 1) {
      message("\nYou may set `nb.cores` to a number higher than 1 for faster computation.\n")
    }

    if (Sys.info()[['sysname']] == "Windows" && nb.cores > 1) {
      message("\nThe argument `nb.cores` is ignore on Windows-based infrastructure.\n")
    }

    time <- system.time({

     h_selection <- "LOO-CV"

     if (length(x_ori) > 500) message("Bandwidth selection using LOO-CV... (may take a while)")

     ## define function for performing leave-one-out cross-validation ## TODO: extract code into standalone function
     CV <- function(h) {
       CVi <- mclapply(seq_along(t_ori), function(oob) {
         obj <- calc_tcor(x = x_ori[-oob], y = y_ori[-oob], t = t_ori[-oob], h = h, t.for.pred = t_ori[oob],
                          cor.method = cor.method, kernel = kernel, param_smoother = param_smoother)
         (obj$r - ((x_ori[oob] - obj$x) * (y_ori[oob] - obj$y)) / (obj$sd_x * obj$sd_y))^2
       })
       CVi_num <- as.numeric(CVi)
       failing <- is.na(CVi_num)
       res <- mean(CVi_num[!failing])
       if (verbose) {
         print(paste("h =", round(h, digits = 1L), "    # failing =", sum(failing),"    CV =", res))
       }
       res
     }

     ## estimate h by cross-validation
     h_max <- 8*sqrt(length(x_ori))
     opt <- stats::optimize(CV, interval = c(1, h_max), tol = 1)

     ## if h_max is best, use elbow criterion instead
     h <- opt$minimum
     CV_error <- opt$objective
     print(paste("h selected using LOO-CV =", round(h, digits = 1L)))
     verbose <- FALSE
     CV_bound <- CV(h_max)

     if (CV_bound <= opt$objective) {
      CV_error <- NA # remove stored value since not used in this case
      h_selection <- "elbow criterion"

      if (length(x_ori) > 500) message("Bandwidth selection using elbow criterion... (may take a while)")

      RMSE_Lh0 <- function(h0) {
        d <- data.frame(h = 1:h0)
        d$CV_h <- sapply(h, CV)
        fit_Lh0 <- stats::lm(CV_h ~ h, data = d)
        sqrt(mean(stats::residuals(fit_Lh0)^2))
      }

      RMSE_Rh0 <- function(h0) {
        d <- data.frame(h = (h0 + 1):h_max)
        d$CV_h <- sapply(h, CV)
        fit_Rh0 <- stats::lm(CV_h ~ h, data = d)
        sqrt(mean(stats::residuals(fit_Rh0)^2))
      }

      RMSE_h0 <- function(h0) {
        (h0 - 1)/(h_max - 1)*RMSE_Lh0(h0) + (h_max - h0)/(h_max - 1)*RMSE_Rh0(h0)
      }

      opt <- stats::optimize(RMSE_h0, interval = c(3, h_max - 1), tol = 1)
      h <- opt$minimum
      print(paste("h selected using elbow criterion =", round(h, digits = 1L)))
     }
   })
  } else {
    time <- NULL
    h_selection <- "fixed by user"
  }

  ## compute correlation with the selected bandwidth
  res_all <- calc_tcor(x = x_ori, y = y_ori, t = t_ori, h = h, cor.method = cor.method,
                       kernel = kernel, param_smoother = param_smoother)

  ## format data for output
  res <- res_all[, c("t", "r")]

  ## restore missing values
  if (keep.missing) {
    res <- merge(res, data.frame(t = t), all.y = TRUE)
  }

  ## store other useful info as attributes
  attr(res, "h") <- h
  attr(res, "CV_error") <- CV_error
  attr(res, "h_selection") <- h_selection
  attr(res, "h_select_duration") <- time
  res
}


#' @describeIn tcor Internal function computing the correlation for a given bandwidth.
#'
#' The function calls the kernel smoothing procedure on each component required
#' to compute the time-varrying correlation. It returns a dataframe with the time,
#' the correlation value and the underlying (smoothed) components.
#'
#' @export
#'
#' @examples
#'
#'
#' ##################################################################
#' ## Examples for the internal function computing the correlation ##
#' ##################################################################
#'
#' with(head(stockprice), calc_tcor(x = SP500, y = FTSE100, t = DateID, h = 20))
#' with(head(stockprice), calc_tcor(x = SP500, y = FTSE100, t = DateID, h = 20,
#'      t.for.pred = DateID[1]))
#'
#' ## The function can handle non consecutive time points
#'
#' set.seed(1)
#' calc_tcor(x = rnorm(10), y = rnorm(10), t = c(1:2, 23:30), h = 2)
#'
#'
#' ## The function can handle non-ordered time series
#'
#' with(head(stockprice)[c(1, 3, 6, 2, 4, 5), ], calc_tcor(x = SP500, y = FTSE100, t = DateID, h = 20))
#'
#'
#' ## Note: the function does not handle missing data (by design)
#'
#' # calc_tcor(x = c(NA, rnorm(9)), y = rnorm(10), t = c(1:2, 23:30), h = 2) ## should err if ran!
#'
calc_tcor <- function(x, y, t = seq_along(x), t.for.pred = t, h, cor.method = c("pearson", "spearman"),
                      kernel = c("epanechnikov", "box", "normal"), param_smoother = list()) {

  ## checking inputs
  if (any(is.na(c(x, y, t)))) {
    stop("Computing correlation requires x, y & t not to contain any NAs; otherwise, r becomes meaningless as it would no longer be bounded between [-1, 1].")
  }

  if (length(x) != length(y)) stop("x and y must have same length")
  if (length(x) != length(t)) stop("t must have same length as x and y")
  #if (length(x) == 0) stop("missing data for requested time(s)")

  if (length(h) != 1) stop("h must be a scalar (numerical value of length 1)")

  cor.method <- match.arg(cor.method)

  ## a spearman correlation equates a pearson on ranked data
  if (cor.method == "spearman") {
    y <- rank(y)
    x <- rank(x)
  }

  ## we create a matrix with the required components
  U <- cbind(x = x, y = y, x2 = x^2, y2 = y^2, xy = x*y)

  ## we smooth each component (i.e. each column in the matrix)
  x_smoothed <- kern_smooth(x = x, t = t, h = h, t.for.pred = t.for.pred,
                            kernel = kernel, param_smoother = param_smoother) ## separate call to retrieve t once (and x)
  other_smoothed_list <- apply(U[, -1], 2L, function(v) {
    kern_smooth(x = v, t = t, h = h, t.for.pred = t.for.pred, kernel = kernel, param_smoother = param_smoother)$x
    }, simplify = FALSE) ## combine call for everything else
  other_smoothed <- do.call("cbind", other_smoothed_list)
  smoothed <- cbind(x_smoothed, as.data.frame(other_smoothed))

  ## we compute the time varrying correlation coefficient
  smoothed$sd_x <- sqrt(smoothed$x2 - smoothed$x^2)
  smoothed$sd_y <- sqrt(smoothed$y2 - smoothed$y^2)
  smoothed$r <- (smoothed$xy - smoothed$x * smoothed$y) / (smoothed$sd_x * smoothed$sd_y)
  smoothed
  }


#' @describeIn tcor Internal function computing the Ht array.
#'
#' Ht is a component needed to compute confindence intervals defined in eq. 6 from Choi & Shin, 2021.
#' The function returns a 5 x 5 x t array.
#'
#' @export
#' @param smoothed_obj an object created with [`calc_tcor`].
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' ####################################################################################
#' ## Examples for the internal function computing Ht (eq. 6 from Choi & Shin, 2021) ##
#' ####################################################################################
#'
#' foo <- with(na.omit(stockprice),
#'             calc_tcor(x = SP500, y = FTSE100, t = DateID, h = 20, kernel = "box"))
#' head(foo)
#' calc_Ht(foo)
#'
calc_Ht <- function(smoothed_obj) {
  res <- array(0, dim = c(5, 5, nrow(smoothed_obj)))

  res[1, 1, ] <- 2*smoothed_obj$x*smoothed_obj$sd_x
  res[1, 2, ] <- 2*smoothed_obj$y*smoothed_obj$sd_y*smoothed_obj$r
  res[1, 3, ] <- smoothed_obj$sd_x
  res[1, 4, ] <- smoothed_obj$sd_y*smoothed_obj$r
  res[1, 5, ] <- smoothed_obj$x*smoothed_obj$sd_y*smoothed_obj$r + smoothed_obj$y*smoothed_obj$sd_x

  res[2, 2, ] <- 2*smoothed_obj$y*smoothed_obj$sd_y*sqrt(1 - smoothed_obj$r^2)
  res[2, 4, ] <- smoothed_obj$sd_y*sqrt(1 - smoothed_obj$r^2)
  res[2, 5, ] <- smoothed_obj$x*smoothed_obj$sd_y*sqrt(1 - smoothed_obj$r^2)

  res[3, 1, ] <- smoothed_obj$sd_x^2
  res[3, 2, ] <- smoothed_obj$sd_y^2*smoothed_obj$r^2
  res[3, 5, ] <- smoothed_obj$sd_x*smoothed_obj$sd_y*smoothed_obj$r

  res[4, 2, ] <- smoothed_obj$sd_y^2*(smoothed_obj$r^2 - 1)

  res[5, 2, ] <- 2*smoothed_obj$sd_y^2*smoothed_obj$r*sqrt(1 - smoothed_obj$r^2)
  res[5, 5, ] <- sqrt(1 - smoothed_obj$r^2)*smoothed_obj$sd_x*smoothed_obj$sd_y

  res
}
