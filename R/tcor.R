#' Compute time varying correlation coefficients
#'
#' The function `tcor` implements (together with its helper function
#' `calc_rho`) the nonparametric estimation of the time varying correlation
#' coefficient proposed by Choi & Shin, 2021. The general idea is to compute a
#' (Pearson) correlation coefficient (`r(x,y) = (mean(xy) - mean(x)*mean(y)) /
#' (sqrt(mean(x^2)-mean(x)^2) * sqrt(mean(y^2)-mean(y)^2))`), but instead of
#' using the means required for such a computation, each component (i.e., `x`,
#' `y`, `x^2`, `y^2`, `x*y`) is smoothed and the smoothed terms are considered
#' in place the original means. The intensity of the smoothing depends on a
#' unique parameter: the bandwidth (`h`). If `h = Inf`, the method produces the
#' original (i.e., time-invariant) correlation value. The smaller the parameter
#' `h`, the more variation in time is being captured. The parameter `h` can be
#' provided by the user; otherwise it is automatically estimated by the internal
#' helper functions `select_h` and `calc_RMSE` (see **Details**).
#'
#' - **Smoothing**: the smoothing of each component is performed by kernel
#' regression. The default is to use the Epanechnikov kernel following Choi &
#' Shin 2021, but other kernels have also been implemented and can thus
#' alternatively be used (see [`kern_smooth`] for details). The normal kernel
#' seems to sometimes lead to very small bandwidth being selected, but the
#' default can lead to numerical issues. We thus recommend always comparing the
#' results from different kernel methods.
#'
#' - **Bandwidth selection**: when the value used to define the bandwidth (`h`)
#' in `tcor` is set to `NULL` (the default), the internal function `select_h`
#' is used to to select the optimal value for `h`. It is first estimated by
#' leave-one-out cross validation (using internally `calc_RMSE`). If the cross
#' validation error (RMSE) is minimal for the maximal value of `h` considered
#' (`8*sqrt(N)`), rather than taking this as the optimal `h` value, the
#' bandwidth becomes estimated using the so-called elbow criterion. This latter
#' method identifies the value `h` after which the cross validation error
#' decreasing very little. The procedure is detailed in section 2.1 in Choi &
#' Shin, 2021.
#'
#' - **Parallel computation**: if `h` is not provided, an automatic bandwidth
#' selection occurs (see above). For large datasets, this step can be
#' computationally demanding. The current implementation thus relies on
#' [`parallel::mclapply`] and is thus only effective for Linux and MacOS.
#' Relying on parallel processing also implies that you call `options("mc.cores"
#' = XX)` beforehand, replacing `XX` by the relevant number of CPU cores you
#' want to use (see **Examples**). For debugging, do use `options("mc.cores" =
#' 1)`, otherwise you may not be able to see the error messages generated in
#' child nodes.
#'
#' - **Confidence interval**: if `CI` is set to `TRUE`, a confidence interval is
#' calculated as described in Choi & Shin, 2021. This is also necessary for using
#' [`equality_test`] to test differences between correlations at two time points.
#' The computation of the confidence intervals involves multiple internal
#' functions (see [`CI`] for details).
#'
#' @inheritParams kern_smooth
#' @param x a numeric vector.
#' @param y a numeric vector of to be correlated with `x`.
#' @param CI a logical specifying if a confidence interval should be computed or not (default = `FALSE`).
#' @param CI.level a scalar defining the level for `CI` (default = 0.95 for 95% CI).
#' @param cor.method a character string indicating which correlation coefficient
#' is to be computed ("pearson", the default; or "spearman").
#' @param keep.missing a logical specifying if time points associated with missing
#' information should be kept in the output (default = `FALSE` to facilitate plotting).
#' @param verbose a logical specifying if information should be displayed to
#' monitor the progress of the cross validation (default = `FALSE`).
#'
#' @return A dataframe containing:
#' - the time points (`t`).
#' - the estimates of the correlation value (`r`).
#'
#' as well as, if `CI = TRUE`:
#' - the Standard Error (`SE`).
#' - the lower boundary of the confidence intervals (`lwr`).
#' - the upper boundary of the confidence intervals (`upr`).
#'
#' Some metadata are also attached to the dataframe (as attributes):
#' - `h` the bandwidth parameter.
#' - `RMSE` the minimal root mean square error when `h` is selected by cross validation.
#' - `h_selection` the method used to select `h`.
#' - `h_select_duration` the computing time spent to select the bandwidth parameter.
#'
#' @name tcor
#' @rdname tcor
#'
#' @references
#' Choi, JE., Shin, D.W. Nonparametric estimation of time varying correlation coefficient.
#' J. Korean Stat. Soc. 50, 333–353 (2021). \doi{10.1007/s42952-020-00073-6}
#'
#' @seealso [`equality_test`], [`kern_smooth`], [`CI`]
#'
#' @importFrom parallel mclapply
#'
NULL


#' @describeIn tcor **The user-level function to be used**.
#'
#' @order 1
#'
#' @export
#'
#' @examples
#'
#' #####################################################
#' ## Examples for the user-level function to be used ##
#' #####################################################
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
#' ## Infinite bandwidth should match fixed correlation coefficients
#' ## nb: `h = Inf` is not supported by default kernel (`kernel = 'epanechnikov'`)
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
#' ## Computing and plotting the confidence interval
#'
#' res_withCI <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = 200, CI = TRUE))
#' with(res_withCI, {
#'      plot(r ~ t, type = "l", ylab = "Cor", xlab = "Time", las = 1, ylim = c(0, 1))
#'      points(lwr ~ t, type = "l", lty = 2)
#'      points(upr ~ t, type = "l", lty = 2)})
#'
#'
#' ## Automatic selection of the bandwidth using parallel processing and comparison
#' ## of the 3 alternative kernels on full dataset
#' # nb: takes a few minutes to run, so not run by default
#'
#' run <- in_pkgdown() || FALSE ## change to TRUE to run the example
#' if (run) {
#'
#' options("mc.cores" = 2L) ## CPU cores to be used for parallel processing
#'
#' res_hauto_epanech <- with(stockprice,
#'          tcor(x = SP500, y = FTSE100, t = DateID, kernel = "epanechnikov")
#'          )
#'
#' res_hauto_box <- with(stockprice,
#'           tcor(x = SP500, y = FTSE100, t = DateID, kernel = "box")
#'           )
#'
#' res_hauto_norm <- with(stockprice,
#'           tcor(x = SP500, y = FTSE100, t = DateID, kernel = "norm")
#'           )
#'
#' plot(res_hauto_epanech, type = "l", col = "red",
#'      ylab = "Cor", xlab = "Time", las = 1, ylim = c(0, 1))
#' points(res_hauto_box, type = "l", col = "blue")
#' points(res_hauto_norm, type = "l", col = "orange")
#' legend("topright", fill = c("red", "blue", "orange"),
#'        legend = c("epanechnikov", "box", "normal"), bty = "n",
#'        title = "Kernel")
#'
#' }
#'
#'
#' ## Comparison of the 3 alternative kernels under same bandwidth
#' ## nb: it requires to have run the previous example
#'
#' if (run) {
#'
#' res_epanech <- with(stockprice,
#'           tcor(x = SP500, y = FTSE100, t = DateID,
#'           kernel = "epanechnikov", h = attr(res_hauto_epanech, "h"))
#'           )
#'
#' res_box <- with(stockprice,
#'            tcor(x = SP500, y = FTSE100, t = DateID,
#'            kernel = "box", h = attr(res_hauto_epanech, "h"))
#'            )
#'
#' res_norm <- with(stockprice,
#'           tcor(x = SP500, y = FTSE100, t = DateID,
#'           kernel = "norm", h = attr(res_hauto_epanech, "h"))
#'           )
#'
#' plot(res_epanech, type = "l", col = "red", ylab = "Cor", xlab = "Time",
#'      las = 1, ylim = c(0, 1))
#' points(res_box, type = "l", col = "grey")
#' points(res_norm, type = "l", col = "orange")
#' legend("topright", fill = c("red", "grey", "orange"),
#'        legend = c("epanechnikov", "box", "normal"), bty = "n",
#'        title = "Kernel")
#'
#' }
#'
#' ## Automatic selection of the bandwidth using parallel processing with CI
#' # nb: takes a few minutes to run, so not run by default
#'
#' run <- in_pkgdown() || FALSE ## change to TRUE to run the example
#' if (run) {
#'
#' res_hauto_epanechCI <- with(stockprice,
#'           tcor(x = SP500, y = FTSE100, t = DateID, CI = TRUE)
#'           )
#'
#' plot(res_hauto_epanechCI[, c("t", "r")], type = "l", col = "red",
#'      ylab = "Cor", xlab = "Time", las = 1, ylim = c(0, 1))
#' points(res_hauto_epanechCI[, c("t", "lwr")], type = "l", col = "red", lty = 2)
#' points(res_hauto_epanechCI[, c("t", "upr")], type = "l", col = "red", lty = 2)
#'
#' }
#'
#'
#' ## Not all kernels work well in all situations
#' ## nb1: EuStockMarkets is a time-series object provided with R
#' ## nb2: takes a few minutes to run, so not run by default
#'
#' run <- in_pkgdown() || FALSE ## change to TRUE to run the example
#' if (run) {
#'
#' EuStock_norm <- tcor(EuStockMarkets[, "DAX"], EuStockMarkets[, "SMI"], kernel = "normal")
#' plot(EuStock_norm, type = "l", las = 1) ## normal kernel seems to work great with these data
#'
#' }
#'
#' run <- in_pkgdown() || FALSE ## change to TRUE to run the example
#' if (run) {
#'
#' EuStock_epanech <- tcor(EuStockMarkets[, "DAX"], EuStockMarkets[, "SMI"])
#' plot(EuStock_epanech, type = "l", las = 1) ## problem of estimation for first and last years
#'
#' }
#'
#'
tcor <- function(x, y, t = seq_along(x), h = NULL, cor.method = c("pearson", "spearman"),
                 kernel = c("epanechnikov", "box", "normal"), CI = FALSE, CI.level = 0.95,
                 param_smoother = list(), keep.missing = FALSE,
                 verbose = FALSE) {

  ## stripping out missing data
  missing <- is.na(x) | is.na(y) | is.na(t)
  x_ori <- x[!missing] ## ori = original, i.e. not smoothed
  y_ori <- y[!missing]
  t_ori <- t[!missing]

  ## if the bandwidth is not given, we must estimate it
  if (is.null(h)) {
    bandwidth_obj <- select_h(x = x_ori, y = y_ori, t = t_ori, cor.method = cor.method, kernel = kernel, param_smoother = param_smoother, verbose = verbose)
  } else {
    bandwidth_obj <- list(h = h, h_selection = "fixed by user", RMSE = NA, time = NULL)
  }


  ## compute correlation with the selected bandwidth
  res_all <- calc_rho(x = x_ori, y = y_ori, t = t_ori, h = bandwidth_obj$h, cor.method = cor.method,
                       kernel = kernel, param_smoother = param_smoother)

  ## format data for output
  res <- res_all[, c("t", "rho_smoothed")]

  ## fix out of range values if required
  out_of_range <- res$rho_smoothed > 1 | res$rho_smoothed < -1
  if (any(out_of_range, na.rm = TRUE)) {
    res$rho_smoothed[res$rho_smoothed > 1] <- 1
    res$rho_smoothed[res$rho_smoothed < -1] <- -1
    warning(paste(sum(out_of_range, na.rm = TRUE), "out of", length(out_of_range), "correlation values were estimated out of the [-1, 1] range and where thus forced to [-1, 1]. Using another kernel may avoid such problem."))
  }

  ## compute CI if requested
  if (CI) {

    res$SE <- calc_SE(smoothed_obj = res_all, h = bandwidth_obj$h, AR.method = "yule-walker") ## for now fixed AR.method
    res$lwr <- res$rho_smoothed + stats::qnorm((1 - CI.level)/2)*res$SE
    res$upr <- res$rho_smoothed + stats::qnorm((1 + CI.level)/2)*res$SE

    ## fix out of range values if required
    out_of_range <- res$lwr > 1 | res$lwr < -1 | res$upr > 1 | res$upr < -1
    if (any(out_of_range, na.rm = TRUE)) {
      res$lwr[res$lwr > 1] <- 1
      res$lwr[res$lwr < -1] <- -1
      res$upr[res$upr > 1] <- 1
      res$upr[res$upr < -1] <- -1
      warning(paste(sum(out_of_range, na.rm = TRUE), "out of", length(out_of_range), "CI boundary values were estimated out of the [-1, 1] range and where thus forced to [-1, 1] Using another kernel may avoid such problem."))
    }
  }

  ## rename column with correlation
  colnames(res)[colnames(res) == "rho_smoothed"] <- "r"

  ## restore missing values
  if (keep.missing) {
    res <- merge(res, data.frame(t = t), all.y = TRUE)
  }

  ## store other useful info as attributes
  attr(res, "call") <- match.call()
  attr(res, "CI") <- CI
  attr(res, "h") <- bandwidth_obj$h
  attr(res, "RMSE") <- bandwidth_obj$RMSE
  attr(res, "h_selection") <- bandwidth_obj$h_selection
  attr(res, "h_select_duration") <- bandwidth_obj$time
  res
}


#' @describeIn tcor Internal function computing the correlation for a given bandwidth.
#'
#' The function calls the kernel smoothing procedure on each component required
#' to compute the time-varying correlation. It returns a dataframe with the time,
#' the correlation value and the underlying components used for the computation.
#'
#' @order 2
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
#' ## Computing the correlation and its component for the first six time points
#'
#' with(head(stockprice), calc_rho(x = SP500, y = FTSE100, t = DateID, h = 20))
#'
#'
#' ## Predicting the correlation and its component at a specific time point
#'
#' with(head(stockprice), calc_rho(x = SP500, y = FTSE100, t = DateID, h = 20,
#'      t.for.pred = DateID[1]))
#'
#'
#' ## The function can handle non consecutive time points
#'
#' set.seed(1)
#' calc_rho(x = rnorm(10), y = rnorm(10), t = c(1:5, 26:30), h = 3, kernel = "box")
#'
#'
#' ## The function can handle non-ordered time series
#'
#' with(head(stockprice)[c(1, 3, 6, 2, 4, 5), ], calc_rho(x = SP500, y = FTSE100, t = DateID, h = 20))
#'
#'
#' ## Note: the function does not handle missing data (by design)
#'
#' # calc_rho(x = c(NA, rnorm(9)), y = rnorm(10), t = c(1:2, 23:30), h = 2) ## should err (if ran)
#'
calc_rho <- function(x, y, t = seq_along(x), t.for.pred = t, h, cor.method = c("pearson", "spearman"),
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
                            kernel = kernel, param_smoother = param_smoother, output = "list") ## separate call to retrieve t once (and x)
  other_smoothed_list <- apply(U[, -1], 2L, function(v) {
    kern_smooth(x = v, t = t, h = h, t.for.pred = t.for.pred, kernel = kernel, param_smoother = param_smoother, output = "list")$x
  }, simplify = FALSE) ## combine call for everything else

  smoothed <- c(x_smoothed, other_smoothed_list) ## output as list, in any case: not to be coerced into matrix -> t can be a Date

  ## we compute the time varying correlation coefficient
  smoothed$sd_x <- sqrt(smoothed$x2 - smoothed$x^2)
  smoothed$sd_y <- sqrt(smoothed$y2 - smoothed$y^2)
  smoothed$rho <- (smoothed$xy - smoothed$x * smoothed$y) / (smoothed$sd_x * smoothed$sd_y)

  ## we rename the components so it is clear they are smoothed
  names(smoothed) <- c(names(smoothed)[1], paste0(names(smoothed)[-1], "_smoothed"))

  ## we add original (non-smoothed) components:
  U_at_t_for_pred <- U[match(t.for.pred, t), , drop = FALSE]
  cbind(U_at_t_for_pred, as.data.frame(smoothed))

}


#' @describeIn tcor Internal function computing the root mean square error (RMSE) for a given bandwidth.
#'
#' The function removes each time point one by one and predicts the correlation
#' at the missing time point based on the other time points. It then computes
#' and returns the RMSE between this predicted correlation and the one predicted
#' using the full dataset. See also *Bandwidth selection* and *Parallel
#' computation* in
#' **Details**.
#'
#' @order 3
#'
#' @export
#'
#' @examples
#'
#'
#' ###########################################################
#' ## Examples for the internal function computing the RMSE ##
#' ###########################################################
#'
#' ## Compute the RMSE on the correlation estimate
#' # nb: takes a few seconds to run, so not run by default
#'
#' run <- in_pkgdown() || FALSE ## change to TRUE to run the example
#' if (run) {
#'
#' small_clean_dataset <- head(na.omit(stockprice), n = 200)
#' with(small_clean_dataset, calc_RMSE(x = SP500, y = FTSE100, t = DateID, h = 10))
#'
#' }
#'
#'
calc_RMSE <- function(h, x, y, t = seq_along(x), cor.method = c("pearson", "spearman"),
                      kernel = c("epanechnikov", "box", "normal"), param_smoother = list(),
                      verbose = FALSE) {

  CVi <- mclapply(seq_along(t), function(oob) { ## TODO: check scheduling effect and try to make this Windows friendly
    obj <- calc_rho(x = x[-oob], y = y[-oob], t = t[-oob], h = h, t.for.pred = t[oob],
                    cor.method = cor.method, kernel = kernel, param_smoother = param_smoother)
    (obj$rho_smoothed - ((x[oob] - obj$x_smoothed) * (y[oob] - obj$y_smoothed)) / (obj$sd_x_smoothed * obj$sd_y_smoothed))^2
  })
  CVi_num <- as.numeric(CVi)
  failing <- is.na(CVi_num)
  res <- sqrt(mean(CVi_num[!failing]))
  if (verbose) {
    print(paste("h =", round(h, digits = 1L), "    # failing =", sum(failing),"    RMSE =", res))
  }

  res
}


#' @describeIn tcor Internal function selecting the optimal bandwidth parameter `h`.
#'
#' The function selects and returns the optimal bandwidth parameter `h` using an
#' optimizer ([stats::optimize()]) which searches the `h` value associated with
#' the smallest RMSE returned by [calc_RMSE()]. See also *Bandwidth selection*
#' in **Details**.
#'
#' @order 4
#'
#' @export
#'
#' @examples
#'
#'
#' ################################################################
#' ## Examples for the internal function selecting the bandwidth ##
#' ################################################################
#'
#' ## Automatic selection of the bandwidth using parallel processing
#' # nb: takes a few seconds to run, so not run by default
#'
#' run <- in_pkgdown() || FALSE ## change to TRUE to run the example
#' if (run) {
#'
#' small_clean_dataset <- head(na.omit(stockprice), n = 200)
#' with(small_clean_dataset, select_h(x = SP500, y = FTSE100, t = DateID))
#'
#' }
#'
select_h <- function(x, y, t = seq_along(x), cor.method = c("pearson", "spearman"),
                     kernel = c("epanechnikov", "box", "normal"), param_smoother = list(),
                     verbose = FALSE) {

  RMSE <- NA # default value for export if not computed

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

  h_selection <- "LOO-CV"

  if (length(x) > 500) message("Bandwidth selection using cross validation... (may take a while)")

  time <- system.time({
    ## estimate h by cross-validation
    h_max <- 8*sqrt(length(x))
    opt <- stats::optimize(calc_RMSE, interval = c(1, h_max),
                           x = x, y = y, t = t, cor.method = cor.method,
                           kernel = kernel, param_smoother = param_smoother,
                           verbose = verbose,
                           tol = 1) ## TODO: check if other optimizer would be best

    ## if h_max is best, use elbow criterion instead
    h <- opt$minimum
    RMSE <- opt$objective
    message(paste("h selected using LOO-CV =", round(h, digits = 1L)))

    RMSE_bound <- calc_RMSE(h = h_max, x = x, y = y, t = t, cor.method = cor.method,
                            kernel = kernel, param_smoother = param_smoother,
                            verbose = FALSE)

    if (RMSE_bound <= opt$objective) {
      RMSE <- NA # remove stored value since not used in this case
      h_selection <- "elbow criterion"

      if (length(x) > 500) message("Bandwidth selection using elbow criterion... (may take a while)")

      calc_RMSE_Lh0 <- function(h0) {
        d <- data.frame(h = 1:h0)
        d$RMSE_h <- sapply(h, calc_RMSE,
                           x = x, y = y, t = t, cor.method = cor.method,
                           kernel = kernel, param_smoother = param_smoother,
                           verbose = FALSE)
        fit_Lh0 <- stats::lm(RMSE_h ~ h, data = d)
        sqrt(mean(stats::residuals(fit_Lh0)^2))
      }

      calc_RMSE_Rh0 <- function(h0) {
        d <- data.frame(h = (h0 + 1):h_max)
        d$RMSE_h <- sapply(h, calc_RMSE,
                           x = x, y = y, t = t, cor.method = cor.method,
                           kernel = kernel, param_smoother = param_smoother,
                           verbose = FALSE)
        fit_Rh0 <- stats::lm(RMSE_h ~ h, data = d)
        sqrt(mean(stats::residuals(fit_Rh0)^2))
      }

      calc_RMSE_h0 <- function(h0) {
        (h0 - 1)/(h_max - 1)*calc_RMSE_Lh0(h0) + (h_max - h0)/(h_max - 1)*calc_RMSE_Rh0(h0)
      }

      opt <- stats::optimize(calc_RMSE_h0, interval = c(3, h_max - 1), tol = 1)
      h <- opt$minimum
      message(paste("h selected using elbow criterion =", round(h, digits = 1L)))
    }
  })
  message(paste("Bandwidth automatic selection completed in", round(time[3], digits = 1L), "seconds"))

  list(h = h, h_selection = h_selection, RMSE = RMSE, time = time[3])
}
