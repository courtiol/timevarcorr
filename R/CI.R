#' Internal functions for the computation of confidence intervals
#'
#' These functions compute the different terms required for [`tcor()`] to compute the confidence
#' interval around the time-varying correlation coefficient. These terms are defined in Choi & Shin (2021).
#'
#' @seealso [`tcor()`]
#' @name CI
#'
#' @references
#' Choi, JE., Shin, D.W. Nonparametric estimation of time varying correlation coefficient.
#' J. Korean Stat. Soc. 50, 333–353 (2021). \doi{10.1007/s42952-020-00073-6}
#'
#' Andrews, D. W. K. Heteroskedasticity and autocorrelation consistent covariance matrix estimation.
#' Econometrica: Journal of the Econometric Society, 817-858 (1991).
#'
#' @return
#' - `calc_H()` returns a 5 x 5 x \eqn{t} array of elements of class numeric, which corresponds to \eqn{\hat{H_t}} in Choi & Shin (2021).
#' - `calc_e()` returns a \eqn{t} x 5 matrix of elements of class numeric storing the residuals, which corresponds to \eqn{\hat{e}_t} in Choi & Shin (2021).
#' - `calc_Gamma()` returns a 5 x 5 matrix of elements of class numeric, which corresponds to \eqn{\hat{\Gamma}_l} in Choi & Shin (2021).
#' - `calc_GammaINF()` returns a 5 x 5 matrix of elements of class numeric, which corresponds to \eqn{\hat{\Gamma}^\infty} in Choi & Shin (2021).
#' - `calc_L_And()` returns a scalar of class numeric, which corresponds to \eqn{L_{And}} in Choi & Shin (2021).
#' - `calc_D()` returns a \eqn{t} x 5 matrix of elements of class numeric storing the residuals, which corresponds to \eqn{D_t} in Choi & Shin (2021).
#' - `calc_SE()` returns a vector of length \eqn{t} of elements of class numeric, which corresponds to \eqn{se(\hat{\rho}_t(h))} in Choi & Shin (2021).
#'
#' @examples
#' rho_obj <- with(na.omit(stockprice),
#'                 calc_rho(x = SP500, y = FTSE100, t = DateID, h = 20, kernel = "box"))
#' head(rho_obj)
#'
NULL


#' @describeIn CI computes the \eqn{\hat{H_t}} array.
#'
#' \eqn{\hat{H_t}} is a component needed to compute confidence intervals;
#' \eqn{H_t} is defined in eq. 6 from Choi & Shin (2021).
#'
#' @export
#' @param smoothed_obj an object created with [`calc_rho`].
#'
#' @examples
#' ## Computing \eqn{\hat{H_t}}
#'
#' H <- calc_H(smoothed_obj = rho_obj)
#' H[, , 1:2] # H array for the first two time points
#'
calc_H <- function(smoothed_obj) {
  res <- array(0, dim = c(5, 5, nrow(smoothed_obj)))

  res[1, 1, ] <- 2*smoothed_obj$x_smoothed*smoothed_obj$sd_x_smoothed
  res[1, 2, ] <- 2*smoothed_obj$y_smoothed*smoothed_obj$sd_y*smoothed_obj$rho_smoothed
  res[1, 3, ] <- smoothed_obj$sd_x_smoothed
  res[1, 4, ] <- smoothed_obj$sd_y_smoothed*smoothed_obj$rho_smoothed
  res[1, 5, ] <- smoothed_obj$x_smoothed*smoothed_obj$sd_y_smoothed*smoothed_obj$rho_smoothed + smoothed_obj$y_smoothed*smoothed_obj$sd_x_smoothed

  res[2, 2, ] <- 2*smoothed_obj$y_smoothed*smoothed_obj$sd_y_smoothed*sqrt(1 - smoothed_obj$rho_smoothed^2)
  res[2, 4, ] <- smoothed_obj$sd_y_smoothed*sqrt(1 - smoothed_obj$rho_smoothed^2)
  res[2, 5, ] <- smoothed_obj$x_smoothed*smoothed_obj$sd_y_smoothed*sqrt(1 - smoothed_obj$rho_smoothed^2)

  res[3, 1, ] <- smoothed_obj$sd_x_smoothed^2
  res[3, 2, ] <- smoothed_obj$sd_y_smoothed^2*smoothed_obj$rho_smoothed^2
  res[3, 5, ] <- smoothed_obj$sd_x_smoothed*smoothed_obj$sd_y_smoothed*smoothed_obj$rho_smoothed

  res[4, 2, ] <- smoothed_obj$sd_y_smoothed^2*(smoothed_obj$rho_smoothed^2 - 1)

  res[5, 2, ] <- 2*smoothed_obj$sd_y_smoothed^2*smoothed_obj$rho_smoothed*sqrt(1 - smoothed_obj$rho_smoothed^2)
  res[5, 5, ] <- sqrt(1 - smoothed_obj$rho_smoothed^2)*smoothed_obj$sd_x_smoothed*smoothed_obj$sd_y_smoothed

  res
}


#' @describeIn CI computes \eqn{\hat{e}_t}.
#'
#' \eqn{\hat{e}_t} is defined in eq. 9 from Choi & Shin (2021).
#'
#' @export
#' @param H an object created with `calc_H`.
#'
#' @examples
#' ## Computing \eqn{\hat{e}_t}
#'
#' e <- calc_e(smoothed_obj = rho_obj, H = H)
#' head(e) # e matrix for the first six time points
#'
calc_e <- function(smoothed_obj, H) {
  res <- matrix(0, ncol = 5, nrow = dim(H)[3])
  for (i in seq_len(dim(H)[3])) {
    Ut   <- t(smoothed_obj[i, c("x2", "y2", "x", "y", "xy"), drop = FALSE])
    mu_t <- t(smoothed_obj[i, c("x2_smoothed", "y2_smoothed", "x_smoothed", "y_smoothed", "xy_smoothed"), drop = FALSE])
    res[i, ] <- tryCatch({solve(t(H[, , i])) %*% (Ut - mu_t)}, error = function(e) stop("The bandwidth parameter `h` is too small. CI cannot be computed."))
  }
  colnames(res) <- paste0(c("x2", "y2", "x", "y", "xy"), "_resid")
  res
}


#' @describeIn CI computes \eqn{\hat{\Gamma}_l}.
#'
#' \eqn{\hat{\Gamma}_l} is defined in eq. 9 from Choi & Shin (2021).
#'
#' @export
#' @param e an object created with `calc_e`.
#' @param l a scalar indicating a number of time points.
#'
#' @examples
#' ## Computing \eqn{\hat{\Gamma}_l}
#'
#' calc_Gamma(e = e, l = 3)
#'
calc_Gamma <- function(e, l) {
  res <- matrix(0, nrow = 5, ncol = 5)
  for (t in seq_len(nrow(e) - l)) {
    res <- res + t(e[t, , drop = FALSE]) %*% e[t + l, , drop = FALSE]
  }
  res/nrow(e)
}


#' @describeIn CI computes \eqn{\hat{\Gamma}^\infty}.
#'
#' \eqn{\hat{\Gamma}^\infty} is the long run variance estimator, defined in eq. 9 from Choi & Shin (2021).
#'
#' @export
#' @param L a scalar indicating a bandwidth parameter.
#'
#' @examples
#' ## Computing \eqn{\hat{\Gamma}^\infty}
#'
#' calc_GammaINF(e = e, L = 2)
#'
calc_GammaINF <- function(e, L) {

  Gamma_0 <- calc_Gamma(e = e, l = 0)

  sum_term <- matrix(0, nrow = 5, ncol = 5)
  for (l in seq_len(L)) {
    Gamma_l <- calc_Gamma(e = e, l = l)
    sum_term <- sum_term + (1 - l/L) * (Gamma_l + t(Gamma_l))
  }

  Gamma_0 + sum_term
}


#' @describeIn CI computes \eqn{L_{And}}.
#'
#' \eqn{L_{And}} is defined in Choi & Shin (2021, p 342).
#' It also corresponds to \eqn{S_T^*}, eq 5.3 in Andrews (1991).
#'
#' @export
#' @param AR.method character string specifying the method to fit the autoregressive model used to compute \eqn{\hat{\gamma}_1} in \eqn{L_{And}} (see [`stats::ar`] for details).
#'
#' @examples
#' ## Computing \eqn{L_{And}}
#'
#' calc_L_And(e = e)
#' sapply(c("yule-walker", "burg", "ols", "mle", "yw"),
#'        function(m) calc_L_And(e = e, AR.method = m)) ## comparing AR.methods
#'
calc_L_And <- function(e, AR.method = c("yule-walker", "burg", "ols", "mle", "yw")) {
  AR.method <- match.arg(AR.method)
  gamma_1 <- stats::ar(rowSums(e), method = AR.method, order.max = 1L)$ar
  if (length(gamma_1) == 0) return(0)
  1.1447*((4*nrow(e)*gamma_1^2)/((1 - gamma_1^2)^2))^(1/3)
}


#' @describeIn CI computes \eqn{D_t}.
#'
#' \eqn{D_t} is defined in Choi & Shin (2021, p 338).
#'
#' @export
#'
#' @examples
#' ## Computing \eqn{D_t}
#'
#' D <- calc_D(smoothed_obj = rho_obj)
#' head(D) # D matrix for the first six time points
#'
calc_D <- function(smoothed_obj) {
  res <- matrix(0, nrow = nrow(smoothed_obj), ncol = 5)
  res[, 1] <- -1/2*(smoothed_obj$xy_smoothed - smoothed_obj$x_smoothed*smoothed_obj$y_smoothed)/(smoothed_obj$sd_x_smoothed^3*smoothed_obj$sd_y_smoothed)
  res[, 2] <- -1/2*(smoothed_obj$xy_smoothed - smoothed_obj$x_smoothed*smoothed_obj$y_smoothed)/(smoothed_obj$sd_x_smoothed*smoothed_obj$sd_y_smoothed^3)
  res[, 3] <- (-smoothed_obj$y_smoothed*smoothed_obj$sd_x_smoothed^2 + smoothed_obj$x_smoothed*(smoothed_obj$xy_smoothed - smoothed_obj$x_smoothed*smoothed_obj$y_smoothed))/(smoothed_obj$sd_x_smoothed^3*smoothed_obj$sd_y_smoothed)
  res[, 4] <- (-smoothed_obj$x_smoothed*smoothed_obj$sd_y_smoothed^2 + smoothed_obj$y_smoothed*(smoothed_obj$xy_smoothed - smoothed_obj$x_smoothed*smoothed_obj$y_smoothed))/(smoothed_obj$sd_x_smoothed*smoothed_obj$sd_y_smoothed^3)
  res[, 5] <- 1/(smoothed_obj$sd_x_smoothed*smoothed_obj$sd_y_smoothed)
  res
}


#' @describeIn CI computes \eqn{se(\hat{\rho}_t(h))}.
#'
#' The standard deviation of the time-varying correlation (\eqn{se(\hat{\rho}_t(h))}) is defined in eq. 8 from Choi & Shin (2021).
#' It depends on \eqn{D_{Lt}}, \eqn{D_{Mt}} & \eqn{D_{Ut}}, themselves defined in Choi & Shin (2021, p 337 & 339).
#' The \eqn{D_{Xt}} terms are all computed within the function since they all rely on the same components.
#'
#' @export
#' @inheritParams kern_smooth
#'
#' @examples
#' ## Computing \eqn{se(\hat{\rho}_t(h))}
#' # nb: takes a few seconds to run
#'
#' run <- FALSE ## change to TRUE to run the example
#' if (in_pkgdown() || run) {
#'
#' SE <- calc_SE(smoothed_obj = rho_obj, h = 50)
#' head(SE) # SE vector for the first six time points
#'
#' }
#'
#'
calc_SE <- function(smoothed_obj, h, AR.method = c("yule-walker", "burg", "ols", "mle", "yw")) {

  ## Bartlett kernel
  bk <- function(x) pmax(0, 1 - abs(x)) ## 1 - abs(x) if abs(x) < 1, 0 otherwise
  bk2 <- function(x) bk(x)^2

  ## Common terms
  H <- calc_H(smoothed_obj = smoothed_obj)
  e <- calc_e(smoothed_obj = smoothed_obj, H = H)
  D <- calc_D(smoothed_obj = smoothed_obj)
  L_And <- calc_L_And(e = e, AR.method = AR.method)
  GammaINF <- calc_GammaINF(e = e, L = L_And)
  N <- nrow(smoothed_obj)

  SE <- numeric(N)
  for (t in seq_len(N)) {
    Dt <- D[t, ]
    Ht <- H[, , t]

    ## boundary cases at the beginning of the time series
    if (t < h) {
      c <- t/h
      DLt <- stats::integrate(bk, -c, 1)$value^-2*t(Ht) %*% GammaINF %*% Ht * stats::integrate(bk2, -c, 1)$value
      SE[t] <- sqrt(t(Dt) %*% DLt %*% Dt/h)
    ## boundary cases at the end of the time series
    } else if (t > N - h) {
      c <- (N - t)/h
      DUt <- stats::integrate(bk, -1, c)$value^-2*t(Ht) %*% GammaINF %*% Ht * stats::integrate(bk2, -1, c)$value
      SE[t] <- sqrt(t(Dt) %*% DUt %*% Dt/h)
    ## non-boundary cases
    } else {
      DMt <- t(Ht) %*% GammaINF %*% Ht * stats::integrate(bk2, -1, 1)$value
      SE[t] <- sqrt(t(Dt) %*% DMt %*% Dt/h)
    }
  }
  SE
}
