#' Internal functions for the computation of confidence intervals
#'
#' These functions compute the different terms required to compute the confidence
#' interval around the time-varying correlation coefficients.
#'
#' @seealso [`tcor`]
#' @name CI
#'
#' @references
#' Choi, JE., Shin, D.W. Nonparametric estimation of time varying correlation coefficient.
#' J. Korean Stat. Soc. 50, 333–353 (2021). https://doi.org/10.1007/s42952-020-00073-6
#'
#' Andrews, D. W. K. Heteroskedasticity and autocorrelation consistent covariance matrix estimation.
#' Econometrica: Journal of the Econometric Society, 817-858 (1991).
#'
#' @examples
#'
#' rho_obj <- with(na.omit(stockprice),
#'                 calc_rho(x = SP500, y = FTSE100, t = DateID, h = 20, kernel = "box"))
#' head(rho_obj)
#'
NULL


#' @describeIn CI Internal function computing the `$\hat{H_t}$` array.
#'
#' `$\hat{H_t}$` is a component needed to compute confidence intervals;
#' `$H_t$` is defined in eq. 6 from Choi & Shin, 2021.
#' The function returns a 5 x 5 x t array.
#'
#' @export
#' @param smoothed_obj an object created with [`calc_rho`].
#'
#' @examples
#' ## Computing `$\hat{H_t}$`
#'
#' H <- calc_H(smoothed_obj = rho_obj)
#' H
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


#' @describeIn CI Internal function computing `$\hat{e}_t$`.
#'
#' `$\hat{e}_t$` is a component needed to compute confidence intervals;
#' it is defined in eq. 9 from Choi & Shin, 2021.
#' The function returns a 5 x t matrix storing the residuals.
#'
#' @export
#' @param H an object created with [`calc_H`].
#'
#' @examples
#' ## Computing `$\hat{e}_t$`
#'
#' e <- calc_e(smoothed_obj = rho_obj, H = H)
#' e
#'
calc_e <- function(smoothed_obj, H) {
  res <- matrix(0, ncol = 5, nrow = dim(H)[3])
  for (i in seq_len(dim(H)[3])) {
    Ut   <- t(smoothed_obj[i, c("x2", "y2", "x", "y", "xy"), drop = FALSE])
    mu_t <- t(smoothed_obj[i, c("x2_smoothed", "y2_smoothed", "x_smoothed", "y_smoothed", "xy_smoothed"), drop = FALSE])
    res[i, ] <- solve(t(H[, , i])) %*% (Ut - mu_t)
  }
  colnames(res) <- paste0(c("x2", "y2", "x", "y", "xy"), "_resid")
  res
}


#' @describeIn CI Internal function computing `$\hat{\Gamma}_l$`.
#'
#' `$\hat{\Gamma}_l$` is a component needed to compute confidence intervals;
#' it is defined in eq. 9 from Choi & Shin, 2021.
#' The function returns a 5 x 5 matrix.
#'
#' @export
#' @param e an object created with [`calc_e`].
#' @param l a scalar indicating a number of time steps.
#'
#' @examples
#' ## Computing `$\hat{\Gamma}_l$`
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


#' @describeIn CI Internal function computing `$\hat{\Gamma}^\Inf$`.
#'
#' `$\hat{\Gamma}^\Inf$` is a component needed to compute confidence intervals (the long run variance estimator);
#' it is defined in eq. 9 from Choi & Shin, 2021.
#' The function returns a 5 x 5 matrix.
#'
#' @export
#' @param L a scalar indicating a bandwidth parameter.
#'
#' @examples
#' ## Computing `$\hat{\Gamma}^\Inf$`
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


#' @describeIn CI Internal function computing `$L_{And}$`.
#'
#' `$L_{And}$` is a component needed to compute confidence intervals;
#' it is defined in Choi & Shin, 2021, p 342.
#' It also corresponds to `$S_T^*$`, eq 5.3 in Andrews 1991.
#' The function returns a scalar which should be used as an input for `L` in [`calc_GammaINF`].
#'
#' @export
#' @param AR.method character string specifying the method to fit the autoregressive model used to compute `$\hat{\gamma}_1$` in `$L_{And}$` (see [`stats::ar`] for details).
#'
#' @examples
#' ## Computing `$L_{And}$`
#'
#' calc_L_And(e = e)
#'
calc_L_And <- function(e, AR.method = c("yule-walker", "burg", "ols", "mle", "yw")) {
  AR.method <- match.arg(AR.method)
  gamma_1 <- stats::ar(rowSums(e), method = AR.method, order.max = 1L)$ar
  1.1447*((4*nrow(e)*gamma_1^2)/((1 - gamma_1^2)^2))^(1/3)
}
